import os, sys, glob
sys.path.append('/home/zhouyj/Xiaojiang/PpkAssocLoc')
import argparse
import numpy as np
import matplotlib.pyplot as plt
from obspy.core import *
from scipy.signal import correlate
import data_pipeline as dp
import pickers
import config
import time
import warnings
warnings.filterwarnings("ignore")

def preprocess(st, decim_rate, freq_band):
    st = st.decimate(decim_rate)
    st = st.detrend('demean')
    st = st.detrend('linear')
    if type(freq_band)==list:
        return st.filter('bandpass', freqmin=freq_band[0], freqmax=freq_band[1])
    elif type(freq_band)==float:
        return st.filter('highpass', freq=freq_band)

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
    norm_data = np.sqrt(data_cum[ntemp:] - data_cum[0:-ntemp])
    norm_data  = norm_data[0:len(cc)-1]#TODO
    if len(norm_data)!=len(cc[1:]): return [0]
    cc = cc[1:]/norm_temp/norm_data
    cc[np.isinf(cc)] = 0.
    cc[np.isnan(cc)] = 0.
    return cc


def get_temp_events():
    f=open(args.temp_pha); lines=f.readlines(); f.close()
    temp_events=[]
    for line in lines:
        info = line.split(',')
        if len(info)==4:
            lat = float(info[1])
            lon = float(info[2])
            dep = float(info[3])
            head= [info[0], lat, lon, dep]
            temp_events.append([head, {}])
        else:
            sta = info[0]
            tp = UTCDateTime(info[1])
            ts = UTCDateTime(info[2])
            temp_events[-1][1][sta] = [tp, ts]
    return temp_events

def plot(cc):
    plt.subplot(4,1,1)
    plt.plot(cc[0])
    plt.subplot(4,1,2)
    plt.plot(cc[1])
    plt.subplot(4,1,3)
    plt.plot(cc[2])
    plt.subplot(4,1,4)
    plt.plot((cc[0]+cc[1]+cc[2])/3.)
    plt.show()


def main(args):

  # i/o file
  if os.path.exists(args.out_ctlg):
      os.unlink(args.out_ctlg)
      os.unlink(args.out_pha)
  out_ctlg = open(args.out_ctlg, 'a')
  out_pha  = open(args.out_pha,  'a')
  temp_events = get_temp_events()
  temp_dir = args.temp_dir
  # MFT params
  cfg = config.Config()
  decim_rate = cfg.decim_rate
  freq_band = cfg.freq_band
  samp_rate = 100./decim_rate
  dt_trig = cfg.dt_trig
  dt_p = cfg.dt_p
  dt_s = cfg.dt_s
  trig_thres = cfg.trig_thres
  mask_len = int(samp_rate *cfg.mask_len)
  picker = pickers.Trad_PS()

  # get time range
  start_date = UTCDateTime(args.time_range.split(',')[0])
  end_date   = UTCDateTime(args.time_range.split(',')[1])
  print('Making augmented catalog')
  print('time range: {} to {}'.format(start_date, end_date))

  # for all days
  num_day = (end_date.date - start_date.date).days
  for day_idx in range(num_day):

    # get data paths
    datetime = start_date + day_idx*86400
    data_dict = dp.get_xj(args.data_dir, datetime)
    if data_dict=={}: continue
    print('-'*40)
    print('detecting %s'%datetime.date)
    print('read data')
    todel=[]
    for sta in data_dict:
        if len(data_dict[sta])!=3:
            todel.append(sta); continue
        st = read(data_dict[sta][0])
        st+= read(data_dict[sta][1])
        st+= read(data_dict[sta][2])
        st = preprocess(st, decim_rate, freq_band)
        data_dict[sta] = st
    for sta in todel: data_dict.pop(sta)

    # run MFT with all templates
    for temp_event in temp_events:
        t0=time.time()
        head = temp_event[0]
        print('-'*40)
        print('template ', head)
        temp_name = head[0]
        ot = UTCDateTime(head[0])
        cc=[]
        picks=[]
        if len([sta for sta in temp_event[1] if sta in data_dict])<4: continue
        for sta, [tp, ts] in temp_event[1].items():
            # get continuous data
            if sta not in data_dict: continue
            st = data_dict[sta]
            # get template data
            temp_paths = os.path.join(temp_dir, temp_name, '*.%s.*'%sta)
            temp_paths = sorted(glob.glob(temp_paths))
            temp = read(temp_paths[0])
            temp+= read(temp_paths[1])
            temp+= read(temp_paths[2])
            temp = preprocess(temp, decim_rate, freq_band)
            # slice template for trigger, P pick and S pick
            temp_trig = temp.slice(tp-dt_trig[0], tp+dt_trig[1])
            temp_p = temp.slice(tp-dt_p[0], tp+dt_p[1])
            temp_s = temp.slice(ts-dt_s[0], ts+dt_s[1])
            # calc trig cc
            cci=[]
            cci.append(calc_cc(st[0].data, temp_trig[0].data))
            cci.append(calc_cc(st[1].data, temp_trig[1].data))
            cci.append(calc_cc(st[2].data, temp_trig[2].data))
            cci_len = min([len(ccij) for ccij in cci])
            cci = [ccij[0:cci_len] for ccij in cci]
            cci = np.sum(cci,axis=0) /3.
            # time shift to ot
            dt_idx = int(samp_rate *(tp-ot-dt_trig[0]))
            if dt_idx>0: cci = cci[dt_idx:]
            else:        cci = np.concatenate((np.zeros(-dt_idx), cci))
            trig_idxs = np.where(cci>trig_thres)[0]
            slide_idx = 0
            num=0
            for _ in trig_idxs:
                num+=1
                # mask cci with max cc
                trig_idx = trig_idxs[trig_idxs>=slide_idx][0]
                cc_max = np.amax(cci[trig_idx:trig_idx+mask_len])
                idx_max= np.argmax(cci[trig_idx:trig_idx+mask_len]) +trig_idx
                oti = datetime + idx_max/samp_rate
                mask_idx0 = max(0, idx_max-mask_len//2)
                mask_idx1 = idx_max+mask_len//2
                cci[mask_idx0:mask_idx1] = cc_max
                # pick P and S by cross-correlation
                tp0 = datetime + (idx_max + dt_idx)/samp_rate
                ts0 = tp0+(ts-tp)
                p_rng = [dt_p[0]+mask_len/samp_rate, dt_p[1]+mask_len/samp_rate]
                s_rng = [dt_s[0]+mask_len/samp_rate, dt_s[1]+mask_len/samp_rate]
                st_p = st.slice(tp0-p_rng[0], tp0+p_rng[1])
                st_s = st.slice(ts0-s_rng[0], ts0+s_rng[1])
                if len(st_p)!=3 or len(st_s)!=3: break
                cc_p  = calc_cc(st_p[2].data, temp_p[2].data)
                cc_s0 = calc_cc(st_s[0].data, temp_s[0].data)
                cc_s1 = calc_cc(st_s[1].data, temp_s[1].data)
                tpi = tp0 -p_rng[0] +dt_p[0] +np.argmax(cc_p)/samp_rate
                tsi = ts0 -s_rng[0] +dt_s[0] +(np.argmax(cc_s0)+np.argmax(cc_s1))/samp_rate/2.
                s_amp0 = picker.get_amp(st_s[0].data)
                s_amp1 = picker.get_amp(st_s[1].data)
                s_amp = (s_amp0+s_amp1)/2.
                picks.append((sta, oti, tpi, tsi, s_amp, np.amax(cc_p), (np.amax(cc_s0)+np.amax(cc_s1))/2.))
                # next trig
                slide_idx = trig_idx + mask_len
                if slide_idx>trig_idxs[-1]: break
            print('{} process {} trigs | time {:.2f}'.format(sta, num, time.time()-t0))
            cc.append(cci)
        # detect with stacked cc trace
        if len(cc)<4: continue
        cc_len = min([len(cci) for cci in cc])
        cc = [cci[0:cc_len] for cci in cc]
        cc_stack = np.sum(cc,axis=0) /len(cc)
        # write catalog
        picks = np.array(picks, dtype=[('sta','O'),
                                       ('ot','O'),
                                       ('tp','O'),
                                       ('ts','O'),
                                       ('s_amp','O'),
                                       ('cc_p','O'),
                                       ('cc_s','O')])
        det_idxs = np.where(cc_stack>trig_thres)[0]
        slide_idx = 0
        num=0
        for _ in det_idxs:
            num+=1
            det_idx = det_idxs[det_idxs>=slide_idx][0]
            if det_idx+2*mask_len>len(cc_stack)-1: break
            idx_max = np.argmax(cc_stack[det_idx:det_idx+2*mask_len])
            cc_max  = np.amax(cc_stack[det_idx:det_idx+2*mask_len])
            oti = datetime + (det_idx+idx_max)/samp_rate
            print('detection: ', oti, round(cc_max,2))
            out_ctlg.write('{},{},{},{},{:.3f}\n'.format(oti, head[1], head[2], head[3], cc_max))
            # write phase
            event_pick = [pick for pick in picks if pick['ot']>oti-mask_len/samp_rate\
                                                and pick['ot']<oti+mask_len/samp_rate]
            out_pha.write('{},{},{},{},{:.3f}\n'.format(oti, head[1], head[2], head[3], cc_max))
            for pick in event_pick:
                out_pha.write('{},{},{},{:.2f},{:.3f},{:.3f}\n'.\
                       format(pick['sta'], pick['tp'], pick['ts'], pick['s_amp'], pick['cc_p'], pick['cc_s']))
            slide_idx = det_idx + 2*mask_len
            if slide_idx>det_idxs[-1]: break
        print('time consumption: {:.2f}'.format(time.time()-t0))
        print('{} detections, {} stations'.format(num, len(temp_event[1])))

  out_ctlg.close()
  out_pha.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--data_dir', type=str,
                        default='/data3/XJ_SAC/[Y-Z]*/*')
    parser.add_argument('--time_range', type=str,
                        default='20160903,20180901')
    parser.add_argument('--temp_dir', type=str,
                        default='./output/ZSY/Templates')
    parser.add_argument('--temp_pha', type=str,
                        default='./output/ZSY/pha_ZSY_temp.dat') 
    parser.add_argument('--out_ctlg', type=str,
                        default='./output/ZSY/catalog_aug_ZSY.dat')
    parser.add_argument('--out_pha', type=str,
                        default='./output/ZSY/phase_aug_ZSY.dat')
    args = parser.parse_args()
    main(args)
