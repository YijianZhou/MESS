import os, sys, glob
sys.path.append('/home/zhouyj/Xiaojiang/PpkDet')
import argparse
import time
import numpy as np
from obspy import read, UTCDateTime
# import functions from PpkDet
import data_pipeline as dp
import pickers
# import MFT functions
from mft_lib import *
import config
import warnings
warnings.filterwarnings("ignore")
import matplotlib.pyplot as plt

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--data_dir', type=str,
                        default='/data3/XJ_SAC/[Y-Z]*/*')
    parser.add_argument('--time_range', type=str,
                        default='20160901,20190201')
    parser.add_argument('--temp_root', type=str,
                        default='./output/zsy/Templates')
    parser.add_argument('--temp_pha', type=str,
                        default='./output/zsy/temp_zsy.pha')
    parser.add_argument('--out_ctlg', type=str,
                        default='./output/zsy/aug_zsy.ctlg')
    parser.add_argument('--out_pha', type=str,
                        default='./output/zsy/aug_zsy.pha')
    args = parser.parse_args()


def preprocess(st):

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


# make temp dict: location + sta obs
def get_temp_dict():

  f=open(args.temp_pha); lines=f.readlines(); f.close()
  temp_dict = {}
  for line in lines:
    info = line.split(',')
    if len(info)==5:
        temp_name = info[0]
        ot   = UTCDateTime(temp_name)
        lat  = float(info[1])
        lon  = float(info[2])
        dep  = float(info[3])
        temp_loc = [ot, lat, lon, dep]
        temp_dict[temp_name] = [temp_loc, {}]
    else:
        sta = info[0]
        tp = UTCDateTime(info[1])
        ts = UTCDateTime(info[2])
        dt_ot = ot - tp + win_trig[0]
        temp_dict[temp_name][1][sta] = [tp, ts, dt_ot]

  print('read template data')
  t=time.time()
  for i,temp_name in enumerate(temp_dict):
    if i%50==0: print(i,'th template', temp_name, int(time.time()-t),'s')
    temp_lov  = temp_dict[temp_name][0]
    [ot, lat, lon, dep] = temp_loc
    temp_picks = temp_dict[temp_name][1]
    for sta in temp_picks:
        [tp, ts, dt_ot] = temp_dict[temp_name][1][sta]
        temp_paths = os.path.join(args.temp_root, temp_name, '*.%s.*'%sta)
        temp_paths = sorted(glob.glob(temp_paths))
        temp = read(temp_paths[0])
        temp+= read(temp_paths[1])
        temp+= read(temp_paths[2])
        temp = preprocess(temp)
        if len(temp)!=3: continue
        temp_trig = temp.slice(tp - win_trig[0], tp + win_trig[1])
        temp_p    = temp.slice(tp - win_p[0], tp + win_p[1])
        temp_s    = temp.slice(ts - win_s[0], ts + win_s[1])
        temp_dict[temp_name][1][sta] = [temp_trig, temp_p, temp_s, tp, ts, ot, dt_ot]
  return temp_dict


# read continuous raw data and preprocess
def read_data(data_dict):

    print('read continuous data')
    todel=[]
    t=time.time()
    for sta in data_dict:
        print(sta, int(time.time()-t),'s')
        if len(data_dict[sta])!=3:
            todel.append(sta); continue
        stream = read(data_dict[sta][0])
        stream+= read(data_dict[sta][1])
        stream+= read(data_dict[sta][2])
        stream = preprocess(stream)
        if len(stream)!=3:
            todel.append(sta); continue
        dt_st = stream[0].stats.starttime - date
        data_dict[sta] = [stream, date, dt_st]
    for sta in todel: data_dict.pop(sta)
    return data_dict



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
temp_dict = get_temp_dict()

# get time range
start_date = UTCDateTime(args.time_range.split(',')[0])
end_date   = UTCDateTime(args.time_range.split(',')[1])
print('Run MFT (cpu version)')
print('time range: {} to {}'.format(start_date, end_date))

# for all days
num_day = (end_date.date - start_date.date).days
for day_idx in range(num_day):

    # get data paths
    date = start_date + day_idx*86400
    data_dict = dp.get_xj(args.data_dir, date)
    if data_dict=={}: continue
    # read data and preprocess
    data_dict = read_data(data_dict)
    print('-'*40)
    print('detecting %s'%date.date)

    # run MFT with all templates
    for temp_name in temp_dict:

        # init
        t0=time.time()
        temp_loc = temp_dict[temp_name][0]
        temp_picks = temp_dict[temp_name][1]
        print('-'*40)
        print('template ', temp_loc)
        if len([sta for sta in temp_picks if sta in data_dict])<4: continue

        # for each station, calc masked cc trace
        cc = calc_masked_cc(temp_picks, data_dict, trig_thres, mask_len, samp_rate)

        # detect with stacked cc trace
        if len(cc)<4: continue
        cc_stack = np.sum(cc,axis=0) / len(cc)
#        plt.plot(cc_stack); plt.show()
        det_ots = det_cc_stack(cc_stack, trig_thres, mask_len, date, samp_rate)
        print('{} detections | time {:.2f}'.format(len(det_ots), time.time()-t0))
        if len(det_ots)==0: continue

        # ppk by cc
        for [det_oti, det_cci] in det_ots:
            picksi = ppk_cc(det_oti, temp_picks, data_dict, 
                            win_p, win_s, picker, mask_len, samp_rate)
            write_det_ppk(det_oti, det_cci, temp_loc[1:], picksi, out_ctlg, out_pha)
        print('time consumption: {:.2f}'.format(time.time()-t0))

out_ctlg.close()
out_pha.close()
