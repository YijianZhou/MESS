""" Cut template data
  Use original catalog as template events, cut into SSD for fast i/o
  Input
    data_dir: dir of continuous data
    pha_path: phase file of original catalog
    out_root: root dir for output
    temp_root: name for template data (root dir)
    chn_codes: channel codes sequence, seperated by comma
  Output
    temp_root/temp_name/net.sta.chn
    Note: temp_name == ot (yyyymmddhhmmss)
"""

import os, sys, glob, shutil
sys.path.append('/home/zhouyj/software/data_prep')
sys.path.append('/home/zhouyj/software/PAD')
import argparse
import numpy as np
from obspy.core import read, UTCDateTime
import sac
import pickers
import data_pipeline as dp
import config

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--data_dir', type=str,
                        default='/data3/luwf_data/Trace/Linear_Pad/*')
    parser.add_argument('--pha_path', type=str,
                        default='/data3/luwf_data/Trace/JZG_temp/jzg_temp.pha')
    parser.add_argument('--out_root', type=str,
                        default='./output/tmp')
    parser.add_argument('--temp_root', type=str,
                        default='Templates')
    parser.add_argument('--chn_codes', type=str,
                        default='HHE,HHN,HHZ')
    args = parser.parse_args()


# i/o files
data_dir = args.data_dir
pha_path = args.pha_path
f=open(pha_path); lines=f.readlines(); f.close()
out_root = args.out_root
temp_root = os.path.join(out_root, args.temp_root)
# cut params
cfg = config.Config()
t_blank = cfg.t_blank # sec before tp (cut win)
win_len = cfg.win_len # sec
chn_codes = args.chn_codes.split(',')


# make phase dict
pha_dict = {}
for line in lines:
  codes = line.split(',')
  if codes[0][0:5]=='2017-': continue
  if len(codes)==5:
    event_name = codes[0]
    ot = UTCDateTime(codes[0])
    lat, lon, dep, mag = [float(code) for code in codes[1:]]
    header = [ot, lat, lon, dep, mag]
    pha_dict[event_name] = [header, []]
  else:
    net, sta = codes[0:2]
    tp, ts = [UTCDateTime(code) for code in codes[2:4]]
    pha_dict[event_name][1].append([net, sta, tp, ts])


# cut all events data
for i,event_name in enumerate(pha_dict):

    # event info
    [header, picks] = pha_dict[event_name]
    ot, lat, lon, dep, mag = header
    data_dict = dp.get_jz(data_dir, ot)
    event_dir = os.path.join(temp_root, event_name)
    if not os.path.exists(event_dir): os.makedirs(event_dir)
    if i%10==0: print('cutting {}th event'.format(i))

    for pick in picks:
        net, sta, tp, ts = pick
        b = tp - UTCDateTime(ot.date) - t_blank
        data_paths = data_dict[sta]
        out_paths = [os.path.join(event_dir,'%s.%s.%s'%(net,sta,chn)) for chn in chn_codes]
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

