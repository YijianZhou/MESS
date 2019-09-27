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
                        default='/data2/Ridgecrest/*/*')
    parser.add_argument('--pha_path', type=str,
                        default='./output/rc/phase_rc1.dat')
    parser.add_argument('--out_root', type=str,
                        default='./output/rc')
    parser.add_argument('--temp_root', type=str,
                        default='Templates')
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


# make phase dict
pha_dict = {}
for line in lines:
  codes = line.split(',')
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


# for all events
for i,event_name in enumerate(pha_dict):

    # event info
    [header, picks] = pha_dict[event_name]
    ot, lat, lon, dep, mag = header
    data_dict = dp.get_rc(data_dir, ot)
    event_dir = os.path.join(temp_root, event_name)
    if not os.path.exists(event_dir): os.makedirs(event_dir)
    if i%10==0: print('cutting {}th event'.format(i))

    for pick in picks:
        net, sta, tp, ts = pick
        b = tp - UTCDateTime(ot.date) - t_blank
        # cut event
        out_paths = []
        for data_path in data_dict[sta]:
            fname = os.path.split(data_path)[-1]
            chn = fname.split('.')[-2] #TODO
            out_paths.append(os.path.join(event_dir, '%s.%s.%s'%(net,sta,chn)))
        sac.cut(data_dict[sta][0], b, b+win_len, out_paths[0])
        sac.cut(data_dict[sta][1], b, b+win_len, out_paths[1])
        sac.cut(data_dict[sta][2], b, b+win_len, out_paths[2])

        # write head
        t0 = t_blank
        t1 = ts -tp + t_blank
        sac.ch_event(out_paths[0], lon, lat, dep, mag, [t0,t1])
        sac.ch_event(out_paths[1], lon, lat, dep, mag, [t0,t1])
        sac.ch_event(out_paths[2], lon, lat, dep, mag, [t0,t1])

