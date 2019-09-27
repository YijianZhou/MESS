import os, sys, glob
sys.path.append('/home/zhouyj/software/PAD')
import argparse
import time
import numpy as np
from obspy import read, UTCDateTime
# import functions from PAD
import data_pipeline as dp
import pickers
# import MFT functions
from dataset_gpu_torch import read_temp, read_data
from mft_lib_gpu_torch import *
import config
import warnings
warnings.filterwarnings("ignore")
import matplotlib.pyplot as plt
import torch.multiprocessing as mp
import torch
mp.set_sharing_strategy('file_system')
os.environ["CUDA_VISIBLE_DEVICES"] = "0"


if __name__ == '__main__':
  mp.set_start_method('spawn')
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
  samp_rate = 100. / decim_rate
  win_p = [int(samp_rate * win) for win in cfg.win_p]
  win_s = [int(samp_rate * win) for win in cfg.win_s]
  min_sta = cfg.min_sta
  trig_thres = cfg.trig_thres
  mask_len = int(samp_rate * cfg.mask_len)
  picker = pickers.Trad_PS()

  # i/o file
  out_ctlg = open(args.out_ctlg,'w')
  out_pha = open(args.out_pha,'w')
  temp_dict = read_temp(args.temp_pha, args.temp_root)

  # get time range
  start_date, end_date = [UTCDateTime(date) for date in args.time_range.split(',')]
  print('Run MFT (gpu version)')
  print('Time range: {} to {}'.format(start_date, end_date))

  # for all days
  num_day = (end_date.date - start_date.date).days
  for day_idx in range(num_day):

    # read data
    date = start_date + day_idx*86400
    data_dict = dp.get_rc(args.data_dir, date)
    if data_dict=={}: continue
    data_dict = read_data(data_dict)
    print('-'*40)
    print('Detecting %s'%date.date)

    # run MFT with all templates
    for temp_name, [temp_loc, pick_dict] in temp_dict.items():

        # init
        t=time.time()
        torch.cuda.empty_cache()
        print('template ', temp_loc)
        # drop bad sta
        todel = [sta for sta in pick_dict if sta not in data_dict]
        for sta in todel: pick_dict.pop(sta)
        if len(pick_dict)<min_sta: continue

        # for each station, calc masked cc trace
        cc_holder = torch.zeros([len(pick_dict), int(86400*samp_rate)])
        cc_masked = calc_masked_cc(cc_holder, pick_dict, data_dict, trig_thres, mask_len)
        # detect with stacked cc trace
        cc_stack = np.sum(cc_masked,axis=0) / len(cc_masked)
#        plt.plot(cc_stack); plt.show()
        det_ots = det_cc_stack(cc_stack, trig_thres, mask_len)
        print('{} detections | time {:.1f}'.format(len(det_ots), time.time()-t))
        if len(det_ots)==0: continue

        # ppk by cc
        print('pick p&s by corr')
        for [det_ot, det_cc] in det_ots:
            picks = ppk_cc(det_ot, pick_dict, data_dict, win_p, win_s, picker, mask_len)
            det_ot = idx2time(det_ot, samp_rate, date)
            for i in range(len(picks)): 
                picks[i][1:3] = [idx2time(idx, samp_rate, date) for idx in picks[i][1:3]]
            write_det_ppk(det_ot, det_cc, temp_name, temp_loc, picks, out_ctlg, out_pha)
        print('time consumption: {:.2f}'.format(time.time()-t))

  out_ctlg.close()
  out_pha.close()
