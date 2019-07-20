import os, glob, time
import torch
from torch.utils.data import Dataset, DataLoader
from obspy import read, UTCDateTime
import numpy as np
import config

# preprocess params
cfg = config.Config()
decim_rate = cfg.decim_rate
samp_rate = 100./decim_rate
freq_band = cfg.freq_band
win_trig = cfg.win_trig
win_p = cfg.win_p
win_s = cfg.win_s
npts_trig = int(sum(win_trig) * samp_rate) + 1
temp_win = [int(sum(win) * samp_rate) + 1 for win in [win_trig, win_p, win_s]]

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

# read 3 chn stream
def read_stream(stream_paths):
    stream = read(stream_paths[0])
    stream+= read(stream_paths[1])
    stream+= read(stream_paths[2])
    return preprocess(stream)


def np2cuda(data):
    return torch.from_numpy(data).cuda()

def st2np(stream):
    return np.array([stream[0].data, stream[1].data, stream[2].data])

def st2cuda(stream):
    return np2cuda(st2np(stream))


class Templates(Dataset):

  def __init__(self, temp_root, temp_dict):
    """ Dataset for reading Templates
    """
    self.temp_dirs = glob.glob(os.path.join(temp_root,'*'))
    self.temp_dict = temp_dict

  def __getitem__(self, index):
    # temp info
    temp_dir = self.temp_dirs[index]
    temp_name = os.path.split(temp_dir)[-1]
    out_dict = {}
    # read data
    temp_paths = glob.glob(os.path.join(temp_dir,'*.*Z'))
    for temp_path in temp_paths:
        is_bad = False
        net, sta, _ = os.path.split(temp_path)[-1].split('.')
        [tp, ts, _] = self.temp_dict[temp_name][1][sta]
        stream_paths = sorted(glob.glob(os.path.join(temp_dir, '%s.%s.*'%(net,sta))))
        stream = read_stream(stream_paths)
        temp_trig = stream.slice(tp - win_trig[0], tp + win_trig[1])
        temp_p    = stream.slice(tp - win_p[0], tp + win_p[1])
        temp_s    = stream.slice(ts - win_s[0], ts + win_s[1])
        temp = [st2cuda(temp_trig), st2np(temp_p), st2np(temp_s)]
        for i,tempi in enumerate(temp):
            if len(tempi[0]) != temp_win[i]: is_bad = True
        if is_bad: continue
        norm_trig = [torch.sqrt(torch.sum(tri**2)) for tri in temp[0]]
        norm_p    = [np.sqrt(np.sum(tri**2)) for tri in temp[1]]
        norm_s    = [np.sqrt(np.sum(tri**2)) for tri in temp[2]]
        norm_temp = [norm_trig, norm_p, norm_s]
        # output (data, label)
        out_dict[sta] = [temp, norm_temp]
    return temp_name, out_dict

  def __len__(self):
    return len(self.temp_dirs)


""" Read Template Data
  Inputs
    temp_pha: phase file (text) for template phases
    temp_root: root dir fot template data
  Outputs
    temp_dict: temp_dict[temp_name] = [temp_loc, temp_picks]
               temp_picks[sta] = [temp, norm_temp] + [tp, ts, ot, dt_ot]
    temp == [temp_trig, temp_p, temp_s]
    norm_temp = [norm_trig, norm_p, norm_s]
"""

def get_temp_dict(temp_pha, temp_root):

    # make temp dict
    f=open(temp_pha); lines=f.readlines(); f.close()
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

    # read temp data
    print('Read template data')
    t=time.time()
    todel = []
    temp_dataset = Templates(temp_root, temp_dict)
    temp_loader = DataLoader(temp_dataset, num_workers=10)
    for i, (temp_name, out_dict) in enumerate(temp_loader):
        temp_name = temp_name[0]
        if len(out_dict)<4: todel.append(temp_name)
        [ot, lat, lon, dep] = temp_dict[temp_name][0]
        for sta in out_dict:
            [tp, ts, dt_ot] = temp_dict[temp_name][1][sta]
            out_dict[sta] += [tp, ts, ot, dt_ot]
        temp_dict[temp_name][1] = out_dict
        if i%100==0: print('{}th template: {} {:.1f}s'.format(i, temp_name, time.time()-t))
    for temp_name in todel: temp_dict.pop(temp_name)
    return temp_dict


""" Read Continuous Data
  Inputs
    data_dict: data_dict[sta] = stream_paths
"""

def get_sta_data(data_paths, date):
    stream = read_stream(data_paths)
    dt_st = stream[0].stats.starttime - date
    # get stream data (np.array)
    st_holder = np.zeros([3,int(86400*samp_rate+1)])
    st_np = st2np(stream.slice(date, date + 86400))
    idx0 = int(dt_st * samp_rate)
    idx1 = int(idx0 + st_np.shape[1])
    st_holder[:, idx0:idx1] = st_np
    st_np = st_holder
    st_cuda = torch.unsqueeze(np2cuda(st_np),0)
    # calc norm data (for corr)
    data_cum = [np.cumsum(datai**2) for datai in st_np]
    norm_data = [np.sqrt(cumi[npts_trig:] - cumi[:-npts_trig]) for cumi in data_cum]
    norm_data = [torch.unsqueeze(np2cuda(normi),0) for normi in norm_data]
    return [[stream, st_cuda], norm_data, date, dt_st]


# read continuous raw data and preprocess
def read_data(data_dict, date):
    t=time.time()
    print('read continuous data')
    num_sta = len(data_dict)
    for sta in data_dict:
        data_dict[sta] = get_sta_data(data_dict[sta], date)
        print('Read {} | time: {:.1f}s'.format(sta, time.time()-t))
    return data_dict

