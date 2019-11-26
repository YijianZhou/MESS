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
min_sta = cfg.min_sta
temp_win = [int(sum(win) * samp_rate) + 1 for win in [win_trig, win_p, win_s]]
npts_trig = temp_win[0]


def preprocess(st):
    # time alignment
    start_time = max([tr.stats.starttime for tr in st])
    end_time   = min([tr.stats.endtime   for tr in st])
    st = st.slice(start_time, end_time)
    # signal process
    st = st.decimate(decim_rate)
    st = st.detrend('demean').detrend('linear').taper(max_percentage=0.05, max_length=10.)
    flt_type = freq_band[0]
    freqmin  = freq_band[1]
    if len(freq_band)==2:
        return st.filter(flt_type, freq=freqmin)
    elif len(freq_band)==3:
        freqmax = freq_band[2]
        return st.filter(flt_type, freqmin=freqmin, freqmax=freqmax)


# read 3 chn stream & prep
def read_stream(stream_paths):
    stream  = read(stream_paths[0])
    stream += read(stream_paths[1])
    stream += read(stream_paths[2])
    return preprocess(stream)


# format trans
def np2cuda(data):
    return torch.from_numpy(data).cuda()

def st2np(stream):
    return np.array([tr.data for tr in stream], dtype=np.float64)

def st2cuda(stream):
    return np2cuda(st2np(stream))



class Templates(Dataset):

  def __init__(self, temp_dict):
    """ Dataset for reading Templates
    """
    self.temp_dict = temp_dict
    self.event_list = sorted(list(temp_dict.keys()))

  def __getitem__(self, index):

    # read one template
    temp_name = self.event_list[index]
    temp_dir, temp_loc, pick_dict = self.temp_dict[temp_name]
    ot = temp_loc[0]

    out_dict = {}
    # read data
    for sta, [tp,ts] in pick_dict.items():
        # read temp stream
        is_bad = False
        stream_paths = sorted(glob.glob(os.path.join(temp_dir, '*.%s.*'%sta)))
        stream = read_stream(stream_paths)
        # cut temp
        temp_trig = stream.slice(tp - win_trig[0], tp + win_trig[1])
        temp_p    = stream.slice(tp - win_p[0],    tp + win_p[1])
        temp_s    = stream.slice(ts - win_s[0],    ts + win_s[1])
        if min([len(st) for st in [temp_trig, temp_p, temp_s]])<3: continue
        temp = [st2np(tempi) for tempi in [temp_trig, temp_p, temp_s]]
        for i,tempi in enumerate(temp):
            if len(tempi[0]) != temp_win[i]: is_bad = True
        if is_bad: continue
        norm_trig = [np.sqrt(np.sum(tr**2)) for tr in temp[0]]
        norm_p    = [np.sqrt(np.sum(tr**2)) for tr in temp[1]]
        norm_s    = [np.sqrt(np.sum(tr**2)) for tr in temp[2]]
        norm_temp = [norm_trig, norm_p, norm_s]
        dt_list = [int(dt*samp_rate) for dt in [tp-ot, ts-ot, ot-tp+win_trig[0]]]
        out_dict[sta] = [temp, norm_temp] + dt_list
    return temp_name, out_dict

  def __len__(self):
    return len(self.event_list)



class Data(Dataset):

  def __init__(self, data_dict):
    """ Dataset for reading Continuous Data
    """
    self.data_dict = data_dict
    self.sta_list = sorted(list(data_dict.keys()))

  def __getitem__(self, index):
    # read one station
    sta = self.sta_list[index]
    # read data
    stream_paths = self.data_dict[sta]
    stream = read_stream(stream_paths)
    start_time = stream[0].stats.starttime
    end_time = stream[0].stats.endtime
    date = UTCDateTime((start_time + (end_time - start_time)/2).date)

    # get stream data (np.array)
    data_holder = np.zeros([3,int(86400*samp_rate+1)])
    data_np = st2np(stream.slice(date, date+86400))
    idx0 = int((start_time - date) * samp_rate)
    idx1 = int(idx0 + data_np.shape[1])
    data_holder[:, idx0:idx1] = data_np
    data_np = data_holder

    # calc norm data (for corr)
    data_cum = [np.cumsum(datai**2) for datai in data_np]
    norm_data = [np.sqrt(cumi[npts_trig:] - cumi[:-npts_trig]) for cumi in data_cum]
    return sta, [data_np, norm_data]

  def __len__(self):
    return len(self.sta_list)



""" Read Template Data
  Inputs
    temp_pha: phase file (text) for template phases
    temp_root: root dir fot template data
  Outputs
    temp_dict = {temp_name: [temp_loc, pick_dict]}
    pick_dict[sta] = [temp, norm_temp] + [ttp, tts, dt_ot]
    temp = [temp_trig, temp_p, temp_s]
    norm_temp = [norm_trig, norm_p, norm_s]
"""

def read_temp(temp_pha, temp_root):

    # 1. read temp pha as temp dict
    print('Reading template phase file')
    f=open(temp_pha); lines=f.readlines(); f.close()
    temp_dict = {}
    for line in lines[0:46]:
      codes = line.split(',')
      if len(codes)==5:
        temp_name = codes[0]
        temp_dir = os.path.join(temp_root, temp_name)
        ot = UTCDateTime(codes[0])
        lat, lon, dep = [float(code) for code in codes[1:4]]
        temp_loc = [ot, lat, lon, dep]
        temp_dict[temp_name] = [temp_dir, temp_loc, {}]
      else:
        sta = codes[1]
        tp, ts = [UTCDateTime(code) for code in codes[2:4]]
        temp_dict[temp_name][-1][sta] = [tp, ts]

    # read temp data
    print('Reading template data')
    t=time.time()
    todel = []
    temp_dataset = Templates(temp_dict)
    temp_loader = DataLoader(temp_dataset, num_workers=10, pin_memory=True)
    for i, (temp_name, pick_dict) in enumerate(temp_loader):
        temp_name = temp_name[0]
        if len(pick_dict)<min_sta: todel.append(temp_name)
        temp_loc = temp_dict[temp_name][1]
        temp_dict[temp_name] = [temp_loc, pick_dict]
        if i%100==0: print('{}th template: {} {:.1f}s'.format(i, temp_name, time.time()-t))
    for temp_name in todel: temp_dict.pop(temp_name)
    return temp_dict


""" Read Continuous Data
  Inputs
    data_dict = {sta: stream_paths}
  Outputs
    data_dict = {sta: [data, norm_data]} 
"""

def read_data(data_dict):
    t=time.time()
    print('Read continuous data')
    data_dataset = Data(data_dict)
    data_loader = DataLoader(data_dataset, num_workers=5, pin_memory=True)
    for i, (sta, outi) in enumerate(data_loader):
        data_dict[sta[0]] = outi
        print('read {} | time {:.1f}s'.format(sta[0], time.time()-t))
    return data_dict

