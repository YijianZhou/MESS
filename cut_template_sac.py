""" Cut template waveform with SAC
  Inputs
    data_dir: dir of continuous data
    temp_pha: template phase file
    out_root: root dir for template data
  Outputs
    temp_root/temp_name/net.sta.chn
    Note: temp_name == ot (yyyymmddhhmmss.ss)
"""
import os, sys, glob, shutil
sys.path.append('/home/zhouyj/software/data_prep')
import argparse
import numpy as np
import multiprocessing as mp
from obspy import read, UTCDateTime
import sac
import config
from dataset_gpu import read_ftemp

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--data_dir', type=str,
                        default='/data3/luwf_data/Trace/Linear_Pad')
    parser.add_argument('--temp_pha', type=str,
                        default='/data3/luwf_data/Trace/JZG_temp/jzg_temp.pha')
    parser.add_argument('--out_root', type=str,
                        default='./output/tmp')
    args = parser.parse_args()


# i/o files
if not os.path.exists(args.out_root): os.makedirs(args.out_root)
temp_list = read_ftemp(args.temp_pha)
# cut params
cfg = config.Config()
num_workers = cfg.num_workers
win_len = cfg.win_len # sec
get_data_dict = cfg.get_data_dict

def cut_event(event_idx):
    # get event info
    id_name, event_loc, pick_dict = temp_list[event_idx]
    event_name = id_name.split('_')[1]
    ot, lat, lon, dep, mag = event_loc
    ot = UTCDateTime(ot)
    data_dict = get_data_dict(ot, args.data_dir)
    event_dir = os.path.join(args.out_root, event_name)
    if not os.path.exists(event_dir): os.makedirs(event_dir)
    # cut event
    print('cutting {}'.format(event_name))
    for net_sta, [tp, ts] in pick_dict.items():
        data_paths = data_dict[net_sta]
        chn_codes = [data_path.split('.')[-2] for data_path in data_paths]
        out_paths = [os.path.join(event_dir,'%s.%s'%(net_sta,chn)) for chn in chn_codes]
        # cut event
        b_list = [tp - read(data_path, headonly=True)[0].stats.starttime - win_len[0] \
            for data_path in data_paths]
        sac.cut(data_paths[0], b_list[0], b_list[0]+sum(win_len), out_paths[0])
        sac.cut(data_paths[1], b_list[1], b_list[1]+sum(win_len), out_paths[1])
        sac.cut(data_paths[2], b_list[2], b_list[2]+sum(win_len), out_paths[2])
        # write header
        t0 = win_len[0]
        t1 = ts - tp + win_len[0]
        tn = {'t0':t0, 't1':t1}
        sac.ch_event(out_paths[0], lon, lat, dep, mag, tn)
        sac.ch_event(out_paths[1], lon, lat, dep, mag, tn)
        sac.ch_event(out_paths[2], lon, lat, dep, mag, tn)

# cut all events data
pool = mp.Pool(num_workers)
pool.map_async(cut_event, range(len(temp_list)))
pool.close()
pool.join()

