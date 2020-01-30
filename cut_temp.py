""" Cut template waveform
  Inputs
    data_dir: dir of continuous data
    temp_pha: phase file of original catalog
    out_root: root dir for output
    temp_root: name for template data (root dir)
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
from dataset_gpu import read_pha, dtime2str

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--data_dir', type=str,
                        default='/data3/luwf_data/Trace/Linear_Pad/*')
    parser.add_argument('--temp_pha', type=str,
                        default='/data3/luwf_data/Trace/JZG_temp/jzg_temp.pha')
    parser.add_argument('--out_root', type=str,
                        default='./output/tmp')
    parser.add_argument('--temp_root', type=str,
                        default='Templates')
    args = parser.parse_args()


# i/o files
temp_root = os.path.join(args.out_root, args.temp_root)
if not os.path.exists(temp_root): os.makedirs(temp_root)
pha_list = read_pha(args.temp_pha)
# cut params
cfg = config.Config()
t_blank = cfg.t_blank # sec before tp (cut win)
win_len = cfg.win_len # sec
chn_dict = cfg.chn_dict
get_data_dict = cfg.get_data_dict

def cut_event(event_id):
    # get event info
    [event_loc, pick_dict] = pha_list[event_id]
    ot, lat, lon, dep, mag = event_loc
    data_dict = get_data_dict(args.data_dir, ot)
    event_name = dtime2str(ot)
    event_dir = os.path.join(temp_root, event_name)
    if not os.path.exists(event_dir): os.makedirs(event_dir)

    # cut event
    print('cutting {}'.format(event_name))
    for net_sta, [tp, ts] in pick_dict.items():
        chn_codes = chn_dict[net_sta.split('.')[0]]
        b = tp - UTCDateTime(ot.date) - t_blank
        data_paths = data_dict[net_sta]
        out_paths = [os.path.join(event_dir,'%s.%s'%(net_sta,chn)) for chn in chn_codes]
        # cut event
        sac.cut(data_paths[0], b, b+win_len, out_paths[0])
        sac.cut(data_paths[1], b, b+win_len, out_paths[1])
        sac.cut(data_paths[2], b, b+win_len, out_paths[2])
        # write header
        t0 = t_blank
        t1 = ts -tp + t_blank
        sac.ch_event(out_paths[0], lon, lat, dep, mag, [t0,t1])
        sac.ch_event(out_paths[1], lon, lat, dep, mag, [t0,t1])
        sac.ch_event(out_paths[2], lon, lat, dep, mag, [t0,t1])

# cut all events data
pool = mp.Pool(processes=10)
pool.map_async(cut_event, range(len(pha_list)))
pool.close()
pool.join()

