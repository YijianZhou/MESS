import os, sys, glob
sys.path.append('/home/zhouyj/Xiaojiang/PpkAssocLoc')
import argparse
import time
import numpy as np
from obspy import read, UTCDateTime
from scipy.signal import correlate
# import functions from PpkAssocLoc
import data_pipeline as dp
import pickers
# import MFT functions
from mft_lib import *
import config
import warnings
warnings.filterwarnings("ignore")


def preprocess(st, decim_rate, freq_band):

    # time alignment
    start_time = max([tr.stats.starttime for tr in st])
    end_time   = min([tr.stats.endtime   for tr in st])
    st = st.slice(start_time, end_time)
    # signal process
    st = st.decimate(decim_rate)
    st = st.detrend('demean').detrend('linear').taper(max_percentage=0.05)
    flt_type = freq_band[0]
    freqmin  = freq_band[1]
    if len(freq_band)==2:
        return st.filter(flt_type, freq=freqmin)
    elif len(freq_band)==3:
        freqmax = freq_band[2]
        return st.filter(flt_type, freqmin=freqmin, freqmax=freqmax)


# get list of template events
def get_temp_events():

    f=open(args.temp_pha); lines=f.readlines(); f.close()
    temp_events = []
    for line in lines:
      info = line.split(',')
      if len(info)==5:
        name = info[0]
        ot   = UTCDateTime(name)
        lat  = float(info[1])
        lon  = float(info[2])
        dep  = float(info[3])
        temp_head= [name, ot, lat, lon, dep]
        temp_events.append([temp_head, {}])
      else:
        sta = info[0]
        tp = UTCDateTime(info[1])
        ts = UTCDateTime(info[2])
        temp_events[-1][1][sta] = [tp, ts]
    return temp_events


# read continuous raw data and preprocess
def read_data(data_dict, decim_rate, freq_band):

    print('read continuous data')
    todel=[]
    for sta in data_dict:
        if len(data_dict[sta])!=3:
            todel.append(sta); continue
        st = read(data_dict[sta][0])
        st+= read(data_dict[sta][1])
        st+= read(data_dict[sta][2])
        st = preprocess(st, decim_rate, freq_band)
        if len(st)!=3:
            todel.append(sta); continue
        data_dict[sta] = st
    for sta in todel: data_dict.pop(sta)
    return data_dict


# read template data and preprocess
def read_temp(temp_events, temp_root, decim_rate, freq_band):

  print('read template data')
  temp_dict = {}
  for temp_event in temp_events:
    temp_name = temp_event[0][0]
    temp_picks = temp_event[1]
    temp_dict[temp_name] = {}
    for sta in temp_picks:
        temp_paths = os.path.join(temp_root, temp_name, '*.%s.*'%sta)
        temp_paths = sorted(glob.glob(temp_paths))
        temp = read(temp_paths[0])
        temp+= read(temp_paths[1])
        temp+= read(temp_paths[2])
        temp = preprocess(temp, decim_rate, freq_band)
        if len(temp)!=3: continue
        temp_dict[temp_name][sta] = temp
  return temp_dict



def main(args):

  # MFT params
  cfg = config.Config()
  decim_rate = cfg.decim_rate
  freq_band = cfg.freq_band
  samp_rate = 100./decim_rate
  win_trig = cfg.win_trig
  win_p = cfg.win_p
  win_s = cfg.win_s
  trig_thres = cfg.trig_thres
  mask_len = int(samp_rate *cfg.mask_len)
  picker = pickers.Trad_PS()

  # i/o file
  out_ctlg = open(args.out_ctlg,'w')
  out_pha  = open(args.out_pha, 'w')
  temp_events = get_temp_events()
  temp_root = args.temp_root
  temp_dict = read_temp(temp_events, temp_root, decim_rate, freq_band)

  # get time range
  start_date = UTCDateTime(args.time_range.split(',')[0])
  end_date   = UTCDateTime(args.time_range.split(',')[1])
  print('Run MFT')
  print('time range: {} to {}'.format(start_date, end_date))

  # for all days
  num_day = (end_date.date - start_date.date).days
  for day_idx in range(num_day):

    # get data paths
    date = start_date + day_idx*86400
    data_dict = dp.get_xj(args.data_dir, date)
    if data_dict=={}: continue
    # read data and preprocess
    data_dict = read_data(data_dict, decim_rate, freq_band)
    print('-'*40)
    print('detecting %s'%date.date)

    # run MFT with all templates
    for temp_event in temp_events:

        # init
        t0=time.time()
        temp_head = temp_event[0]
        print('-'*40)
        print('template ', temp_head)
        temp_name = temp_head[0]
        temp_ot = temp_head[1]
        if len([sta for sta in temp_event[1] if sta in data_dict])<4: continue

        # for each station, calc masked cc trace
        cc = calc_masked_cc(temp_event, temp_dict, data_dict, 
                            win_trig, trig_thres, mask_len, date, samp_rate)

        # detect with stacked cc trace
        if len(cc)<4: continue
        cc_stack = np.sum(cc,axis=0) / len(cc)
        det_ots = det_cc_stack(cc_stack, trig_thres, mask_len, date, samp_rate)
        print('{} detections | time {:.2f}'.format(len(det_ots), time.time()-t0))
        if len(det_ots)==0: continue

        # ppk by cc
        for [det_oti, det_cci] in det_ots:
            picksi = ppk_cc(det_oti, temp_event, temp_dict, data_dict, 
                            win_p, win_s, picker, mask_len, samp_rate)
            write_det_ppk(det_oti, det_cci, temp_head[2:], picksi, out_ctlg, out_pha)
        print('time consumption: {:.2f}'.format(time.time()-t0))

  out_ctlg.close()
  out_pha.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--data_dir', type=str,
                        default='/data3/XJ_SAC/*/*')
    parser.add_argument('--time_range', type=str,
                        default='20160901,20190201')
    parser.add_argument('--temp_root', type=str,
                        default='./output/ZSY/Templates')
    parser.add_argument('--temp_pha', type=str,
                        default='./output/ZSY/pha_ZSY_temp.dat') 
    parser.add_argument('--out_ctlg', type=str,
                        default='./output/ZSY/catalog_aug_ZSY.dat')
    parser.add_argument('--out_pha', type=str,
                        default='./output/ZSY/phase_aug_ZSY.dat')
    args = parser.parse_args()
    main(args)
