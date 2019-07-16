import os, sys, glob
sys.path.append('/home/zhouyj/Documents/PpkDet')
import argparse
import time
import numpy as np
from obspy import read, UTCDateTime
# import functions from PpkDet
import data_pipeline as dp
import pickers
# import MFT functions
from dataset import get_temp_dict, read_data
from mft_lib import *
import config
import warnings
warnings.filterwarnings("ignore")
import matplotlib.pyplot as plt

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--data_dir', type=str,
                        default='/data/ZSY_SAC/[Y-Z]*/*')
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
temp_dict = get_temp_dict(args.temp_pha, args.temp_root)

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
    data_dict = read_data(data_dict, date)
    print('-'*40)
    print('detecting %s'%date.date)

    # run MFT with all templates
    for temp_name in temp_dict:

        # init
        t=time.time()
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
        print('{} detections | time {:.2f}'.format(len(det_ots), time.time()-t))
        if len(det_ots)==0: continue

        # ppk by cc
        for [det_oti, det_cci] in det_ots:
            picksi = ppk_cc(det_oti, temp_picks, data_dict, 
                            win_p, win_s, picker, mask_len, samp_rate)
            write_det_ppk(det_oti, det_cci, temp_loc[1:], picksi, out_ctlg, out_pha)
        print('time consumption: {:.2f}'.format(time.time()-t))

out_ctlg.close()
out_pha.close()
