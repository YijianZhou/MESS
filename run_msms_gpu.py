import os, sys, glob
import argparse
import numpy as np
from obspy import UTCDateTime
# import MSMS functions
from dataset_gpu import read_temp, read_data
from msms_lib_gpu import msms_det, corr_ppk, write_ctlg, write_pha
import config
# filter warnings
import warnings
warnings.filterwarnings("ignore")
# torch for GPU
import torch.multiprocessing as mp
import torch
mp.set_sharing_strategy('file_system')
os.environ["CUDA_VISIBLE_DEVICES"] = "0"

if __name__ == '__main__':
  mp.set_start_method('spawn', force=True) # or 'forkserver'
  parser = argparse.ArgumentParser()
  parser.add_argument('--data_dir', type=str,
                      default='/data/*')
  parser.add_argument('--date_range', type=str,
                      default='20170927-20170928')
  parser.add_argument('--temp_root', type=str,
                      default='./output/Templates')
  parser.add_argument('--temp_pha', type=str,
                      default='./output/temp.pha')
  parser.add_argument('--out_ctlg', type=str,
                      default='./output/tmp.ctlg')
  parser.add_argument('--out_pha', type=str,
                      default='./output/tmp.pha')
  args = parser.parse_args()


  # MSMS params
  cfg = config.Config()
  min_sta = cfg.min_sta

  # i/o paths
  out_ctlg = open(args.out_ctlg,'w')
  out_pha = open(args.out_pha,'w')
  temp_dict = read_temp(args.temp_pha, args.temp_root)

  # get date range
  start_date, end_date = [UTCDateTime(date) for date in args.date_range.split('-')]
  print('run MSMS (gpu version)')
  print('date range: {} to {}'.format(start_date, end_date))

  # for all days
  num_day = (end_date.date - start_date.date).days
  for day_idx in range(num_day):
    # read data
    torch.cuda.empty_cache()
    date = start_date + day_idx*86400
    print('-'*40)
    print('detecting %s'%date.date)
    data_dict = read_data(date, args.data_dir)
    if len(data_dict)<min_sta: continue

    # for all templates
    for temp_name, [temp_loc, temp_pick_dict] in temp_dict.items():
        # get pick 
        print('template name (ot):', temp_name)
        todel = [net_sta for net_sta in temp_pick_dict if net_sta not in data_dict]
        for net_sta in todel: temp_pick_dict.pop(net_sta)
        if len(temp_pick_dict)<min_sta: continue

        # msms det
        dets = msms_det(temp_pick_dict, data_dict)
        if len(dets)==0: continue
        # corr ppk
        for [det_ot, det_cc] in dets:
            picks = corr_ppk(det_ot, temp_pick_dict, data_dict)
            det_ot = date + det_ot
            print('det_ot {}, det_cc {:.2f}'.format(det_ot, det_cc))
            for i in range(len(picks)):
                picks[i][1:3] = [date + dt for dt in picks[i][1:3]]
                print('{0[0]} {0[1]} {0[2]} {0[4]:.3f},{0[5]:.3f}'.format(picks[i]))
            write_ctlg(det_ot, det_cc, temp_name, temp_loc, out_ctlg)
            write_pha(det_ot, det_cc, temp_name, temp_loc, picks, out_pha)

  out_ctlg.close()
  out_pha.close()
