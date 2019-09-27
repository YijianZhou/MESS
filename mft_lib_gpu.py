import torch
import torch.nn.functional as F
from scipy.signal import correlate
import numpy as np
import matplotlib.pyplot as plt
import time

def calc_cc_gpu(data_mat, temp_mat, norm_data_mat, norm_temp_vec):
    """ Cala CC trace by GPU cuda
    """
    num_sta, len_data = data_mat.shape
    _, len_temp = temp_mat.shape
    data_mat = data_mat.view([1, num_sta, len_data])
    temp_mat = temp_mat.view([num_sta, 1, len_temp])
    cc = F.conv1d(data_mat, temp_mat, groups=num_sta)[0]
    cc = cc[:,1:] /norm_data_mat /norm_temp_vec.view(num_sta,1)
    cc[torch.isinf(cc)] = 0.
    cc[torch.isnan(cc)] = 0.
    return cc


def calc_cc(data, temp, norm_data=None, norm_temp=None):
    """ Calc CC trace for a template trace
    Input
        data: continuous data to detect from, np.array
        temp: template data, np.array
    Output
        cc: cross-correlation function
    """
    ntemp, ndata = len(temp), len(data)
    if ntemp>ndata: return [0]
    if not norm_temp: 
        norm_temp = np.sqrt(np.sum(temp**2))
    if not norm_data:
        data_cum = np.cumsum(data**2)
        norm_data = np.sqrt(data_cum[ntemp:] - data_cum[:-ntemp])
    cc = correlate(data, temp, mode='valid')
    cc = cc[1:] /norm_temp /norm_data
    cc[np.isinf(cc)] = 0.
    cc[np.isnan(cc)] = 0.
    return cc


def calc_masked_cc(cc_holder, pick_dict, data_dict, trig_thres, mask_len):

    t=time.time()
    # 1. calc cc (GPU)
    # make input matrix
    data_list, temp_list, dt_list = [], [], []
    for sta, [temp, norm_temp, _, _, dt_ot] in pick_dict.items():
        # get data
        data_tensor, norm_data = data_dict[sta]
        data_list.append([data_tensor, norm_data])
        temp_list.append([temp[0], norm_temp[0]])
        dt_list.append(dt_ot)
    cc = calc_shifted_cc(cc_holder, data_list, temp_list, dt_list)
    cc = cc.cpu().numpy()

    # 2. mask cc traces
    cc_masked = [mask_cc(cci, trig_thres, mask_len) for cci in cc]
    print('process {} stations | time {:.1f}s'.format(len(dt_list), time.time()-t))
    return cc_masked


# 1. calc cc trace & time shift
def calc_shifted_cc(cc_holder, data_list, temp_list, dt_list):

  for i in range(3):
    # data_mat
    data_mat = torch.cat([datai[:,i] for [datai,_] in data_list]).cuda()
    temp_mat = torch.cat([tempi[:,i] for [tempi,_] in temp_list]).cuda()
    norm_data = torch.cat([normi[i] for [_,normi] in data_list]).cuda()
    norm_temp = torch.tensor([normi[i] for [_,normi] in temp_list]).cuda()
    # calc cc traces (mat) with GPU
    if i==0: cc  = calc_cc_gpu(data_mat, temp_mat, norm_data, norm_temp)
    else:    cc += calc_cc_gpu(data_mat, temp_mat, norm_data, norm_temp)
  cc /= 3.

  # time shift to ot
  for i,dt_ot in enumerate(dt_list):
    cci = cc[i][max(0,-dt_ot) : cc_holder.shape[1] - dt_ot]
    cc_holder[i][max(0,dt_ot) : max(0,dt_ot) + len(cci)] = cci
  return cc_holder


# 2. mask cc trace
def mask_cc(cci, trig_thres, mask_len):
  # cc mask
  trig_idxs = np.where(cci > trig_thres)[0]
  slide_idx = 0
  for _ in trig_idxs:
    # mask cci with max cc
    trig_idx = trig_idxs[trig_idxs >= slide_idx][0]
    cc_max = np.amax(cci[trig_idx : trig_idx + mask_len])
    idx_max = np.argmax(cci[trig_idx : trig_idx + mask_len])
    idx_max += trig_idx # to abs idx
    mask_idx0 = max(0, idx_max - mask_len //2)
    mask_idx1 = idx_max + mask_len//2
    cci[mask_idx0 : mask_idx1] = cc_max
    # next trig
    slide_idx = idx_max + 2*mask_len
    if slide_idx > trig_idxs[-1]: break
  return cci


# 3. detect in stacked cc trace
def det_cc_stack(cc_stack, trig_thres, mask_len):

  det_idxs = np.where(cc_stack > trig_thres)[0]
  slide_idx = 0
  det_ots = []
  for _ in det_idxs:

    # this detection
    det_idx = det_idxs[det_idxs >= slide_idx][0]
    if det_idx + 2*mask_len > len(cc_stack)-1: break

    # pick ot
    cc_max  = np.amax(cc_stack[det_idx : det_idx + 2*mask_len])
    idx_max = np.where(cc_stack[det_idx: det_idx + 2*mask_len] == cc_max)
    idx_max = int(np.median(idx_max)) + det_idx # to abs idx
    det_ot = idx_max
    det_ots.append([det_ot, cc_max])

    # next detection
    slide_idx = idx_max + 2*mask_len
    if slide_idx > det_idxs[-1]: break
  return det_ots


# 4. pick P&S by cc
def ppk_cc(det_ot, pick_dict, data_dict, win_p, win_s, picker, mask_len):

  picks = []
  for sta, [temp, norm_temp, ttp, tts, _] in pick_dict.items():

    # get data
    st_tensor, _ = data_dict[sta]
    # org tp & ts (in idx)
    tp0, ts0 = int(det_ot+ttp), int(det_ot+tts)
    # cut p&s data (by points)
    p_rng = [win_p[0]+  mask_len, win_p[1]+  mask_len]
    s_rng = [win_s[0]+2*mask_len, win_s[1]+2*mask_len]
    if tp0 < p_rng[0] or ts0 < s_rng[0]\
    or ts0 + s_rng[1] > st_tensor.shape[-1]: continue
    st_p = st_tensor[0,:, tp0-p_rng[0] : tp0+p_rng[1]].numpy()
    st_s = st_tensor[0,:, ts0-s_rng[0] : ts0+s_rng[1]].numpy()

    # ppk by cc
    temp_p, temp_s = temp[1][0].numpy(), temp[2][0].numpy()
    cc_p  = calc_cc(st_p[2], temp_p[2], None, norm_temp[1][2].numpy())
    cc_s0 = calc_cc(st_s[0], temp_s[0], None, norm_temp[2][0].numpy())
    cc_s1 = calc_cc(st_s[1], temp_s[1], None, norm_temp[2][1].numpy())
    # tp & ts (in idx)
    tp = tp0 + np.argmax(cc_p) - mask_len
    ts = ts0 + (np.argmax(cc_s0) + np.argmax(cc_s1))/2. - 2*mask_len
    cc_p_max = np.amax(cc_p)
    cc_s_max = (np.amax(cc_s0) + np.amax(cc_s1)) / 2.

    # get S amplitude
    amp = [picker.get_amp(tr)**2 for tr in st_s]
    s_amp = np.sqrt(sum(amp))
    # add pick
    picks.append([sta, tp, ts, s_amp, cc_p_max, cc_s_max])
  return picks


def write_det_ppk(det_ot, det_cc, temp_name, temp_loc, picks, out_ctlg, out_pha):

    # write catalog
    _, lat, lon, dep = temp_loc
    event_line = '{},{},{},{},{},{:.3f}\n'.\
    format(temp_name, det_ot, lat, lon, dep, det_cc)
    print('detection: ', event_line[:-1])
    out_ctlg.write(event_line)
    # write phase
    out_pha.write(event_line)
    for pick in picks:
        out_pha.write('{0[0]},{0[1]},{0[2]},{0[3]},{0[4]:.3f},{0[5]:.3f}\n'.format(pick))

def idx2time(idx, samp_rate, date):
    return date + idx / samp_rate
