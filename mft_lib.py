from scipy.signal import correlate
import numpy as np

def calc_cc(data, temp):

    """ Calc CC trace for a template trace
    Input
        data: continuous data to detect from, np.array
        temp: template data, np.array
    Output
        cc: cross-correlation function
    """

    ntemp = len(temp)
    ndata = len(data)
    cc = correlate(data, temp, mode='valid')
    norm_temp = np.sqrt(np.sum(temp**2))
    data_cum = np.cumsum(data**2)
    norm_data = np.sqrt(data_cum[ntemp:] - data_cum[:-ntemp])
    cc = cc[1:] /norm_temp /norm_data
    cc[np.isinf(cc)] = 0.
    cc[np.isnan(cc)] = 0.
    return cc


# 1. calc cc trace & time shift to ot
def calc_shifted_cc(temp_trig, st, dt, samp_rate):

    # calc trig cc
    cci = calc_cc(st[0].data, temp_trig[0].data)
    cci+= calc_cc(st[1].data, temp_trig[1].data)
    cci+= calc_cc(st[2].data, temp_trig[2].data)
    cci = cci / 3.

    # time shift to ot
    cci_holder = np.zeros(int(86400*samp_rate))
    dt = int(samp_rate * dt)
    cci = cci[max(0,-dt) : max(0,-dt) + len(cci)]
    cci_holder[max(0,dt) : max(0,dt) + len(cci)] = cci
    return cci_holder


# 2. calc masked cc trace
def calc_masked_cc(cci, trig_thres, mask_len, temp_p, temp_s, st, win_p, win_s, date, dt_ot_tp, dt_s_p, samp_rate, sta, picker):

  # cc mask
  trig_idxs = np.where(cci > trig_thres)[0]
  slide_idx = 0
  picksi = []
  for _ in trig_idxs:

    # mask cci with max cc
    trig_idx = trig_idxs[trig_idxs >= slide_idx][0]
    cc_max = np.amax(cci[trig_idx : trig_idx + mask_len])
    idx_max = np.argmax(cci[trig_idx : trig_idx + mask_len]) + trig_idx
    oti = date + idx_max / samp_rate
    mask_idx0 = max(0, idx_max - mask_len //2)
    mask_idx1 = idx_max + mask_len//2
    cci[mask_idx0 : mask_idx1] = cc_max

    # pick P and S by cross-correlation
    tpi0 = date + (idx_max + dt_ot_tp) /samp_rate
    tsi0 = tpi0 + dt_s_p
    # cut p&s data
    p_rng = [win_p[0] +  mask_len /samp_rate,
             win_p[1] +  mask_len /samp_rate]
    s_rng = [win_s[0] +2*mask_len /samp_rate,
             win_s[1] +2*mask_len /samp_rate]
    st_p = st.slice(tpi0 - p_rng[0], tpi0 + p_rng[1])
    st_s = st.slice(tsi0 - s_rng[0], tsi0 + s_rng[1])
    if len(st_p)!=3 or len(st_s)!=3: break
    # ppk by cc
    cc_p  = calc_cc(st_p[2].data, temp_p[2].data)
    cc_s0 = calc_cc(st_s[0].data, temp_s[0].data)
    cc_s1 = calc_cc(st_s[1].data, temp_s[1].data)
    tpi = tpi0 - p_rng[0] + win_p[0] + np.argmax(cc_p) /samp_rate
    tsi = tsi0 - s_rng[0] + win_s[0] + (np.argmax(cc_s0) + np.argmax(cc_s1)) /samp_rate/2.
    cc_p_max = np.amax(cc_p)
    cc_s_max = (np.amax(cc_s0) + np.amax(cc_s1)) /2

    # get S amplitude
    ampx = picker.get_amp(st_s[0].data)
    ampy = picker.get_amp(st_s[1].data)
    ampz = picker.get_amp(st_s[2].data)
    s_amp = np.sqrt(ampx**2 + ampy**2 + ampz**2)
    picksi.append((sta, oti, tpi, tsi, s_amp, cc_p_max, cc_s_max))
    # next trig
    slide_idx = trig_idx + mask_len
    if slide_idx>trig_idxs[-1]: break
  return cci, picksi


# 3. detect in stacked cc trace
def det_cc_stack(cc_stack, trig_thres, mask_len, date, samp_rate):

  det_idxs = np.where(cc_stack > trig_thres)[0]
  slide_idx = 0
  ots = []
  for _ in det_idxs:

    # this detection
    det_idx = det_idxs[det_idxs >= slide_idx][0]
    if det_idx + 2*mask_len > len(cc_stack)-1: break

    # pick ot
    cc_max  = np.amax(cc_stack[det_idx : det_idx + 2*mask_len])
    idx_max = np.where(cc_stack[det_idx: det_idx + 2*mask_len] == cc_max)
    idx_max = int(np.median(idx_max))
    oti = date + (det_idx + idx_max) /samp_rate
    ots.append([oti, cc_max])
    print('detection: ', oti, round(cc_max,2))

    # next detection
    slide_idx = det_idx + 2*mask_len
    if slide_idx>det_idxs[-1]: break
  return ots


def write_det_ppk(ots, temp_head, picks, ot_dev, out_ctlg, out_pha):

  for [oti, cc_max] in ots:
    # write catalog
    loc = temp_head[1:]
    out_ctlg.write('{},{},{},{},{:.3f}\n'.format(oti,loc[0],loc[1],loc[2],cc_max))
    # write phase
    event_pick = [pick for pick in picks if pick['ot'] > oti - ot_dev\
                                        and pick['ot'] < oti + ot_dev]
    out_pha.write('{},{},{},{},{:.3f}\n'.format(oti,loc[0],loc[1],loc[2],cc_max))
    for pick in event_pick:
        out_pha.write('{},{},{},{:.2f},{:.3f},{:.3f}\n'.\
        format(pick['sta'],pick['tp'],pick['ts'],pick['s_amp'],pick['cc_p'],pick['cc_s']))

