import os, time
import torch
import torch.nn.functional as F
from scipy.signal import correlate
import numpy as np
from numba import jit
from dataset_gpu import cpu2cuda
import config

# MSMS params
cfg = config.Config()
freq_band = cfg.freq_band
samp_rate = cfg.samp_rate
temp_win_p = [int(samp_rate * win) for win in cfg.temp_win_p]
temp_win_s = [int(samp_rate * win) for win in cfg.temp_win_s]
ppk_win_p = [int(samp_rate * win) for win in cfg.ppk_win_p]
ppk_win_s = [int(samp_rate * win) for win in cfg.ppk_win_s]
mask_len = int(samp_rate * cfg.mask_len)
det_gap = int(samp_rate * cfg.det_gap)
trig_thres = cfg.trig_thres
picker = cfg.picker


def msms_det(temp_pick_dict, data_dict):
    """ MSMS detection (main)
    Inputs
      temp, norm_temp, dt_list = temp_pick_dict[net_sta]
      data, norm_data = data_dict[net_sta]
    Outputs
      dets = [det_ot, det_cc]
    """
    t=time.time()
    # prep input
    num_sta = len(temp_pick_dict)
    cc_holder = np.zeros([num_sta, int(86400*samp_rate)])
    data_list, temp_list, dt_ot_list = [], [], []
    for net_sta, [temp, norm_temp, dt_list] in temp_pick_dict.items():
        data, norm_data = data_dict[net_sta][1:3]
        data_list.append([data, norm_data])
        temp_list.append([temp[0], norm_temp[0]])
        dt_ot_list.append(dt_list[0])
    # 1. match
    cc_mat = match_filter(data_list, temp_list)
    # 2. shift
    cc = shift_ot(cc_holder, cc_mat, dt_ot_list)
    # 3. mask
    cc_masked = [mask_cc(cci) for cci in cc]
    # 4. stack & detect
    cc_stack = np.sum(cc_masked, axis=0) / len(cc_masked)
    dets = det_cc_stack(cc_stack)
    print('{} dets, {} sta, {:.1f}s'.format(len(dets), num_sta, time.time()-t))
    return dets


def corr_ppk(det_ot, temp_pick_dict, data_dict):
    """ Cross-correlation phase picking (main)
    Inputs
      det_ot (sec): detected orgin time (relative)
      temp_pick_dict, data_dict
    Outputs
      picks (list): [net_sta, tp, ts, s_amp, cc_p, cc_s]
    """
    det_ot *= samp_rate # sec to idx
    picks = []
    for net_sta, [temp, norm_temp, dt_list] in temp_pick_dict.items():
        # get np data & temp
        data_np = data_dict[net_sta][0].numpy()
        temp = [tempi.numpy() for tempi in temp]
        norm_temp = [normi.numpy() for normi in norm_temp]
        # slice p&s data
        p_rng = [temp_win_p[0] + ppk_win_p[0], temp_win_p[1] + ppk_win_p[1]]
        s_rng = [temp_win_s[0] + ppk_win_s[0], temp_win_s[1] + ppk_win_s[1]]
        tp0, ts0 = int(det_ot + dt_list[1]), int(det_ot + dt_list[2])
        if tp0 < p_rng[0] or ts0 < s_rng[0] \
        or ts0 + s_rng[1] > data_np.shape[-1]: continue
        data_p = data_np[:, tp0 - p_rng[0] : tp0 + p_rng[1]]
        data_s = data_np[:, ts0 - s_rng[0] : ts0 + s_rng[1]]

        # 1. ppk by cc
        cc_p  = calc_cc(data_p[2], temp[1][2], norm_temp=norm_temp[1][2])
        cc_s0 = calc_cc(data_s[0], temp[2][0], norm_temp=norm_temp[2][0])
        cc_s1 = calc_cc(data_s[1], temp[2][1], norm_temp=norm_temp[2][1])
        # tp & ts (rel sec)
        tp = (tp0 + np.argmax(cc_p) - ppk_win_p[0]) / samp_rate
        ts = (ts0 + (np.argmax(cc_s0) + np.argmax(cc_s1))/2. - ppk_win_s[0]) / samp_rate
        # cc_p & cc_s
        cc_p_max = np.amax(cc_p)
        cc_s_max = (np.amax(cc_s0) + np.amax(cc_s1)) / 2.

        # 2. get amplitude
        amp_xyz = np.array([picker.get_amp(tr) for tr in data_s])
        s_amp = np.linalg.norm(amp_xyz)
        picks.append([net_sta, tp, ts, s_amp, cc_p_max, cc_s_max])
    return picks


""" Base functions
"""

def calc_cc_gpu(data_mat, temp_mat, norm_data_mat, norm_temp_vec):
    num_sta, len_data = data_mat.shape
    _,       len_temp = temp_mat.shape
    data_mat = data_mat.view([1, num_sta, len_data])
    temp_mat = temp_mat.view([num_sta, 1, len_temp])
    cc_mat = F.conv1d(data_mat, temp_mat, groups=num_sta)[0,:,1:]
    cc_mat /= norm_data_mat
    cc_mat /= norm_temp_vec.view([num_sta,1])
    return cc_mat


def calc_cc(data, temp, norm_data=None, norm_temp=None):
    ntemp, ndata = len(temp), len(data)
    if ntemp>ndata: return [0]
    if not norm_temp:
        norm_temp = np.sqrt(np.sum(temp**2))
    if not norm_data:
        data_cum = np.cumsum(data**2)
        norm_data = np.sqrt(data_cum[ntemp:] - data_cum[:-ntemp])
    cc = correlate(data, temp, mode='valid')[1:]
    cc /= norm_data * norm_temp
    cc[np.isinf(cc)] = 0.
    cc[np.isnan(cc)] = 0.
    return cc


# 1. matched filter (calc cc traces)
def match_filter(data_list, temp_list):
    for i in range(3):
        data_mat = torch.stack([datai[i] for [datai,_] in data_list])
        norm_data = torch.stack([normi[i] for [_,normi] in data_list])
        temp_mat = cpu2cuda(torch.stack([tempi[i] for [tempi,_] in temp_list]))
        norm_temp = cpu2cuda(torch.tensor([normi[i] for [_,normi] in temp_list]))
        if i==0: cc_mat  = calc_cc_gpu(data_mat, temp_mat, norm_data, norm_temp)
        else:    cc_mat += calc_cc_gpu(data_mat, temp_mat, norm_data, norm_temp)
    cc_mat = (cc_mat/3.).cpu().numpy()
    cc_mat[np.isinf(cc_mat)] = 0.
    cc_mat[np.isnan(cc_mat)] = 0.
    return cc_mat


# 2. shift time shift to ot
@jit
def shift_ot(cc_holder, cc_mat, dt_ot_list):
    for i,dt_ot in enumerate(dt_ot_list):
        cci = cc_mat[i][max(0,-dt_ot) : cc_holder.shape[1] - dt_ot]
        cc_holder[i][max(0,dt_ot) : max(0,dt_ot) + len(cci)] = cci
    return cc_holder


# 3. mask cc trace with peak values
@jit
def mask_cc(cc):
    trig_idxs = np.where(cc>trig_thres)[0]
    slide_idx = 0
    for trig_idx in trig_idxs:
        if trig_idx < slide_idx: continue
        # mask cc with peak ccs
        cc_trig = cc[trig_idx : trig_idx+2*mask_len]
        cc_max = np.amax(cc_trig)
        idx_max = trig_idx + np.argmax(cc_trig)
        idx0 = max(0, idx_max - mask_len//2)
        idx1 = idx_max + mask_len//2
        cc[idx0:idx1] = cc_max
        # next trig
        slide_idx = trig_idx + 2*mask_len + det_gap
    return cc


# 4. detect on stacked cc trace
@jit
def det_cc_stack(cc_stack):
    det_idxs = np.where(cc_stack>trig_thres)[0]
    slide_idx = 0
    dets = []
    for det_idx in det_idxs:
        if det_idx < slide_idx: continue
        # det ot (rel sec)
        cc_det = cc_stack[det_idx : det_idx+2*mask_len]
        cc_max = np.amax(cc_det)
        det_ot = (det_idx + np.median(np.where(cc_det == cc_max)[0])) / samp_rate
        dets.append([det_ot, cc_max]) 
        # next det
        slide_idx = det_idx + 2*mask_len + det_gap
    return dets


# write detection to catalog
def write_ctlg(det_ot, det_cc, temp_name, temp_loc, out_ctlg):
    out_ctlg.write('{0},{1},{2[1]},{2[2]},{2[3]},{3:.3f}\n'\
        .format(temp_name, det_ot, temp_loc, det_cc))

# write phase picks to phase file
def write_pha(det_ot, det_cc, temp_name, temp_loc, picks, out_pha):
    out_pha.write('{0},{1},{2[1]},{2[2]},{2[3]},{3:.3f}\n'\
        .format(temp_name, det_ot, temp_loc, det_cc))
    for pick in picks:
        out_pha.write('{0[0]},{0[1]},{0[2]},{0[3]},{0[4]:.3f},{0[5]:.3f}\n'.format(pick))

