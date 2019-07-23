import torch
import torch.nn.functional as F
from scipy.signal import correlate
import multiprocessing as mp
import numpy as np
import time

def calc_cc_gpu(data_mat, temp_mat, norm_data_mat, norm_temp_vec):

    num_sta, len_data = data_mat.shape
    _, len_temp = temp_mat.shape
    data_mat = data_mat.view([1, num_sta, len_data])
    temp_mat = temp_mat.view([num_sta, 1, len_temp])
    cc = F.conv1d(data_mat, temp_mat, groups=num_sta)[0]
    cc = cc[:,1:] /norm_data_mat /norm_temp_vec.view(num_sta,1)
    cc[torch.isinf(cc)] = 0.
    cc[torch.isnan(cc)] = 0.
    return cc


def calc_cc(data, temp, **kwargs):

    """ Calc CC trace for a template trace
    Input
        data: continuous data to detect from, np.array
        temp: template data, np.array
    Output
        cc: cross-correlation function
    """
    ntemp = len(temp)
    ndata = len(data)
    if ntemp>ndata: return [0]
    if 'norm_temp' in kwargs: norm_temp = kwargs['norm_temp']
    else: norm_temp = np.sqrt(np.sum(temp**2))
    data_cum = np.cumsum(data**2)
    norm_data = np.sqrt(data_cum[ntemp:] - data_cum[:-ntemp])
    cc = correlate(data, temp, mode='valid')
    cc = cc[1:] /norm_temp /norm_data
    cc[np.isinf(cc)] = 0.
    cc[np.isnan(cc)] = 0.
    return cc


def calc_masked_cc(temp_picks, data_dict, trig_thres, mask_len, samp_rate):

    t=time.time()
    # 1. calc cc (GPU)
    # make input matrix
    data_list = []
    temp_list = []
    dt_list = []
    for sta, [temp, norm_temp, _, _, _, dt_ot] in temp_picks.items():
        # get data
        if sta not in data_dict: continue
        [st_tensor, norm_data, dt_st, _] = data_dict[sta]
        data_list.append([st_tensor, norm_data])
        temp_list.append([temp[0], norm_temp[0]])
        dt_list.append(dt_st+dt_ot)
    cc = calc_shifted_cc(data_list, temp_list, dt_list, samp_rate)
    cc = cc.cpu().numpy()
    # 2. mask cc traces
    cc_masked = []
    for cci in cc:
        cci, num = mask_cc(cci, trig_thres, mask_len)
        cc_masked.append(cci)
    print('process {} stations | time {:.1f}s'.format(len(temp_picks), time.time()-t))
    return cc_masked


# 1. calc cc trace & time shift
def calc_shifted_cc(data_list, temp_list, dt_list, samp_rate):

    # data_mat
    st_mat_e = torch.cat([datai[:,0] for [datai,_] in data_list]).cuda()
    st_mat_n = torch.cat([datai[:,1] for [datai,_] in data_list]).cuda()
    st_mat_z = torch.cat([datai[:,2] for [datai,_] in data_list]).cuda()
    # temp_mat
    temp_mat_e = torch.cat([tempi[:,0] for [tempi,_] in temp_list]).cuda()
    temp_mat_n = torch.cat([tempi[:,1] for [tempi,_] in temp_list]).cuda()
    temp_mat_z = torch.cat([tempi[:,2] for [tempi,_] in temp_list]).cuda()
    # norm_data mat
    norm_data_e = torch.cat([normi[0] for [_,normi] in data_list]).cuda()
    norm_data_n = torch.cat([normi[1] for [_,normi] in data_list]).cuda()
    norm_data_z = torch.cat([normi[2] for [_,normi] in data_list]).cuda()
    # norm_temp vec
    norm_temp_e = torch.tensor([normi[0] for [_,normi] in temp_list]).cuda()
    norm_temp_n = torch.tensor([normi[1] for [_,normi] in temp_list]).cuda()
    norm_temp_z = torch.tensor([normi[2] for [_,normi] in temp_list]).cuda()

    # calc cc traces (mat) with GPU
    cc  = calc_cc_gpu(st_mat_e, temp_mat_e, norm_data_e, norm_temp_e)
    cc += calc_cc_gpu(st_mat_n, temp_mat_n, norm_data_n, norm_temp_n)
    cc += calc_cc_gpu(st_mat_z, temp_mat_z, norm_data_z, norm_temp_z)
    cc /= 3.

    # time shift to ot
    cc_out = []
    cc_holder = torch.zeros([len(dt_list), int(86400*samp_rate)])
    for i,dt in enumerate(dt_list):
        cci = cc[i]
        dt = int(samp_rate * dt)
        cci = cci[max(0,-dt) : max(0,-dt) + len(cci)]
        cc_holder[i][max(0,dt) : max(0,dt) + len(cci)] = cci[0 : cc_holder.shape[1] - max(0,dt)]
    return cc_holder


# 2. mask cc trace
def mask_cc(cci, trig_thres, mask_len):

  # cc mask
  trig_idxs = np.where(cci > trig_thres)[0]
  slide_idx = -1
  num=0
  for _ in trig_idxs:
    num+=1

    # mask cci with max cc
    trig_idx = trig_idxs[trig_idxs > slide_idx][0]
    cc_max = np.amax(cci[trig_idx : trig_idx + mask_len])
    idx_max = np.argmax(cci[trig_idx : trig_idx + mask_len]) + trig_idx
    mask_idx0 = max(0, idx_max - mask_len //2)
    mask_idx1 = idx_max + mask_len//2
    cci[mask_idx0 : mask_idx1] = cc_max

    # next trig
    slide_idx = idx_max + 2*mask_len
    if slide_idx >= trig_idxs[-1]: break
  return cci, num


# 3. detect in stacked cc trace
def det_cc_stack(cc_stack, trig_thres, mask_len, date, samp_rate):

  det_idxs = np.where(cc_stack > trig_thres)[0]
  slide_idx = -1
  det_ots = []
  for _ in det_idxs:

    # this detection
    det_idx = det_idxs[det_idxs > slide_idx][0]
    if det_idx + 2*mask_len > len(cc_stack)-1: break

    # pick ot
    cc_max  = np.amax(cc_stack[det_idx : det_idx + 2*mask_len])
    idx_max = np.where(cc_stack[det_idx: det_idx + 2*mask_len] == cc_max)
    idx_max = int(np.median(idx_max)) + det_idx
    det_oti = date + (det_idx + idx_max) / samp_rate
    det_ots.append([det_oti, cc_max])
    print('detection: ', det_oti, round(cc_max,2))

    # next detection
    slide_idx = idx_max + 2*mask_len
    if slide_idx >= det_idxs[-1]: break
  return det_ots


# 4. pick P&S by cc
def ppk_cc(det_oti, temp_picks, data_dict, 
           win_p, win_s, picker, mask_len, samp_rate):

  picksi = []
  for sta, [temp, norm_temp, tp, ts, ot, _] in temp_picks.items():

    # get data
    if sta not in data_dict: continue
    [st_tensor, _, _, date] = data_dict[sta]

    tpi0 = int(samp_rate * (det_oti + (tp-ot) - date))
    tsi0 = int(samp_rate * (det_oti + (ts-ot) - date))

    # cut p&s data (by points)
    p_rng = [int(win_p[0] * samp_rate) +   mask_len,
             int(win_p[1] * samp_rate) +   mask_len]
    s_rng = [int(win_s[0] * samp_rate) + 2*mask_len,
             int(win_s[1] * samp_rate) + 2*mask_len]
    if tpi0<p_rng[0] or tsi0<s_rng[0]\
    or tsi0+s_rng[1] > st_tensor.shape[-1]: continue
    st_p = st_tensor[0,:, tpi0-p_rng[0] : tpi0+p_rng[1]].numpy()
    st_s = st_tensor[0,:, tsi0-s_rng[0] : tsi0+s_rng[1]].numpy()

    # ppk by cc
    temp_p, temp_s = temp[1][0].numpy(), temp[2][0].numpy()
    cc_p  = calc_cc(st_p[2], temp_p[2], norm_temp=norm_temp[1][2].numpy())
    cc_s0 = calc_cc(st_s[0], temp_s[0], norm_temp=norm_temp[2][0].numpy())
    cc_s1 = calc_cc(st_s[1], temp_s[1], norm_temp=norm_temp[2][1].numpy())
    tpi = tpi0 - p_rng[0] + np.argmax(cc_p) + win_p[0]*samp_rate
    tsi = tsi0 - s_rng[0] + (np.argmax(cc_s0) + np.argmax(cc_s1))/2. + win_s[0]*samp_rate
    tpi = date + tpi / samp_rate
    tsi = date + tsi / samp_rate
    cc_p_max = np.amax(cc_p)
    cc_s_max = (np.amax(cc_s0) + np.amax(cc_s1)) / 2.

    # get S amplitude
    ampx = picker.get_amp(st_s[0])
    ampy = picker.get_amp(st_s[1])
    ampz = picker.get_amp(st_s[2])
    s_amp = np.sqrt(ampx**2 + ampy**2 + ampz**2)

    # add pick
    picksi.append([sta, tpi, tsi, s_amp, cc_p_max, cc_s_max])
  return picksi


def write_det_ppk(det_oti, det_cci, temp_loc, picks, out_ctlg, out_pha):

    # write catalog
    event_line = '{},{},{},{},{:.3f}\n'.\
    format(det_oti, temp_loc[0], temp_loc[1], temp_loc[2], det_cci)
    out_ctlg.write(event_line)
    # write phase
    out_pha.write(event_line)
    for pick in picks:
        out_pha.write('{},{},{},{:.2f},{:.3f},{:.3f}\n'.\
        format(pick[0], pick[1], pick[2], pick[3], pick[4], pick[5]))

