""" Cut template waveform with SAC
  Inputs
    data_dir: dir of continuous data
    temp_pha: template phase file
    out_root: root dir for template data
  Outputs
    temp_root/temp_name/net.sta.chn
    Note: temp_name == ot (yyyymmddhhmmss.ss)
"""
import os, glob, shutil
import argparse
import numpy as np
import torch.multiprocessing as mp
from torch.utils.data import Dataset, DataLoader
from obspy import read, UTCDateTime
import config
from dataset_gpu import read_ftemp, preprocess
import subprocess
os.putenv("SAC_DISPLAY_COPYRIGHT", '0')
import warnings
warnings.filterwarnings("ignore")

# cut params
cfg = config.Config()
num_workers = cfg.num_workers
win_len = cfg.win_len 
get_data_dict = cfg.get_data_dict


def sac_cut(fpath, b, e, out_path):
    p = subprocess.Popen(['sac'], stdin=subprocess.PIPE)
    s = "wild echo off \n"
    s += "cuterr fillz \n"
    s += "cut %s %s \n" %(b, e)
    s += "r %s \n" %fpath
    s += "ch allt (0-&1,b&) iztype IB \n"
    s += "w %s \n" %out_path
    s += "q \n"
    p.communicate(s.encode())


class Cut_Templates(Dataset):
  """ Dataset for cutting templates
  """
  def __init__(self, temp_list):
    self.temp_list = temp_list
    self.data_dir = args.data_dir
    self.out_root = args.out_root

  def __getitem__(self, index):
    data_paths_i = []
    # get event info
    id_name, event_loc, pick_dict = self.temp_list[index]
    event_name = id_name.split('_')[1]
    ot, lat, lon, dep, mag = event_loc
    ot = UTCDateTime(ot)
    data_dict = get_data_dict(ot, self.data_dir)
    event_dir = os.path.join(self.out_root, event_name)
    if not os.path.exists(event_dir): os.makedirs(event_dir)
    # cut event
    for net_sta, [tp, ts] in pick_dict.items():
        data_paths = data_dict[net_sta]
        chn_codes = [data_path.split('.')[-2] for data_path in data_paths]
        out_paths = [os.path.join(event_dir,'%s.%s'%(net_sta,chn)) for chn in chn_codes]
        # cut event
        b_list = [tp - read(data_path, headonly=True)[0].stats.starttime - win_len[0] \
            for data_path in data_paths]
        sac_cut(data_paths[0], b_list[0], b_list[0]+sum(win_len), out_paths[0])
        sac_cut(data_paths[1], b_list[1], b_list[1]+sum(win_len), out_paths[1])
        sac_cut(data_paths[2], b_list[2], b_list[2]+sum(win_len), out_paths[2])
        # preprocess
        st  = read(out_paths[0])
        st += read(out_paths[1])
        st += read(out_paths[2])
        st = preprocess(st)
        if len(st)!=3: 
            for out_path in out_paths: os.unlink(out_path)
            continue
        # write header & record out_paths
        t0 = win_len[0]
        t1 = ts - tp + win_len[0]
        for ii in range(3): 
            st[ii].stats.sac.t0, st[ii].stats.sac.t1 = t0, t1
            st[ii].write(out_paths[ii], format='sac')
        data_paths_i.append(out_paths)
    return data_paths_i

  def __len__(self):
    return len(self.temp_list)


if __name__ == '__main__':
    mp.set_start_method('spawn', force=True) # 'spawn' or 'forkserver'
    parser = argparse.ArgumentParser()
    parser.add_argument('--data_dir', type=str,
                        default='/data/Example_data')
    parser.add_argument('--temp_pha', type=str,
                        default='input/example.temp')
    parser.add_argument('--out_root', type=str,
                        default='output/example_templates')
    args = parser.parse_args()

    # i/o files
    if not os.path.exists(args.out_root): os.makedirs(args.out_root)
    temp_list = read_ftemp(args.temp_pha)
    # for sta-date pairs
    data_paths  = []
    dataset = Cut_Templates(temp_list)
    dataloader = DataLoader(dataset, num_workers=num_workers, batch_size=None)
    for i, data_paths_i in enumerate(dataloader):
        data_paths += data_paths_i
        if i%10==0: print('%s/%s events done/total'%(i,len(dataset)))
    fout_data_paths = os.path.join(args.out_root,'data_paths.npy')
    np.save(fout_data_paths, data_paths)

