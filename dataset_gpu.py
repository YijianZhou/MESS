import os, glob, time
import torch
from torch.utils.data import Dataset, DataLoader
from obspy import read, UTCDateTime
import numpy as np
import config

# params
cfg = config.Config()
get_data_dict = cfg.get_data_dict
num_workers = cfg.num_workers
samp_rate = cfg.samp_rate
freq_band = cfg.freq_band
temp_win_trig = cfg.temp_win_trig
temp_win_p = cfg.temp_win_p
temp_win_s = cfg.temp_win_s
min_sta = cfg.min_sta
temp_win_npts = [int(sum(win) * samp_rate) for win in [temp_win_trig, temp_win_p, temp_win_s]]
npts_trig =  temp_win_npts[0] # for data norm


def read_data(date, data_dir):
    """ Read continuous data
    Inputs
      data_dict = {sta: stream_paths}
    Outputs
      data_dict = {sta: [[data_cuda, norm_data_cuda], [data, norm_data]]} 
    """
    t=time.time()
    print('reading continuous data')
    data_dict = get_data_dict(data_dir, date)
    data_dataset = Data(data_dict)
    data_loader = DataLoader(data_dataset, num_workers=num_workers, batch_size=None, pin_memory=True)
    for i, (sta, datai) in enumerate(data_loader):
        data_dict[sta] = [datai[0]] + [cpu2cuda(di) for di in datai]
        print('read {} | time {:.1f}s'.format(sta, time.time()-t))
    return data_dict


def read_temp(temp_pha, temp_root):
    """ Read templates
    Inputs
      temp_pha (txt): template phase file, evt + pha lines
        event line: ot, lat, lon, dep, mag
        phase line: net, sta, tp, ts, s_amp, p_snr, s_snr
      temp_root: root dir for template data
        temp_root/temp_name/net.sta.chn
        *note: temp_name == ot (yyyymmddhhmmss.ss)
    Outputs
      temp_dict = {temp_name: [temp_loc, temp_pick_dict]}
      temp_pick_dict[net_sta] = [temp, norm_temp, dt_list]
      temp = [temp_trig, temp_p, temp_s]
      norm_temp = [norm_trig, norm_p, norm_s]
      dt_list = [ttp, tts, dt_ot]
    """
    # 1. read phase file
    print('reading template phase file')
    f=open(temp_pha); lines=f.readlines(); f.close()
    temp_dict = {}
    for line in lines:
        codes = line.split(',')
        if len(codes)==5:
            temp_name = codes[0]
            temp_dir = os.path.join(temp_root, temp_name)
            ot = UTCDateTime(codes[0])
            lat, lon, dep = [float(code) for code in codes[1:4]]
            temp_loc = [ot, lat, lon, dep]
            temp_dict[temp_name] = [temp_dir, temp_loc, {}]
        else:
            net_sta = '.'.join(codes[0:2])
            tp, ts = [UTCDateTime(code) for code in codes[2:4]]
            temp_dict[temp_name][-1][net_sta] = [tp, ts]

    # 2. read temp data
    print('reading templates')
    t=time.time()
    todel = []
    temp_dataset = Templates(temp_dict)
    temp_loader = DataLoader(temp_dataset, num_workers=num_workers, batch_size=None, pin_memory=True)
    for i, (temp_name, temp_pick_dict) in enumerate(temp_loader):
        if len(temp_pick_dict)<min_sta: todel.append(temp_name)
        temp_loc = temp_dict[temp_name][1]
        temp_dict[temp_name] = [temp_loc, temp_pick_dict]
        if i%100==0: print('{}th template | time {:.1f}s'.format(i, time.time()-t))
    for temp_name in todel: temp_dict.pop(temp_name)
    return temp_dict


class Data(Dataset):
  """ Dataset for reading continuous data
  """
  def __init__(self, data_dict):
    self.data_dict = data_dict
    self.sta_list = sorted(list(data_dict.keys()))

  def __getitem__(self, index):
    # read stream
    sta = self.sta_list[index]
    stream_paths = self.data_dict[sta]
    stream = read_stream(stream_paths)
    start_time = stream[0].stats.starttime
    end_time = stream[0].stats.endtime
    date = UTCDateTime((start_time + (end_time - start_time)/2).date)
    stream = trim_stream(stream, date, date+86400)
    data_np = st2np(stream)[:, 0:86400*samp_rate]

    # calc norm data (for corr)
    data_cum = [np.cumsum(datai**2) for datai in data_np]
    norm_data = np.array([np.sqrt(cumi[npts_trig:] - cumi[:-npts_trig]) for cumi in data_cum])
    return sta, [data_np, norm_data]

  def __len__(self):
    return len(self.sta_list)


class Templates(Dataset):
  """ Dataset for reading templates
  """
  def __init__(self, temp_dict):
    self.temp_dict = temp_dict
    self.temp_list = sorted(list(temp_dict.keys()))

  def __getitem__(self, index):
    # read one template
    temp_name = self.temp_list[index]
    temp_dir, temp_loc, pick_dict_ppk = self.temp_dict[temp_name]
    ot = temp_loc[0]

    # read data
    pick_dict_data = {}
    for net_sta, [tp,ts] in pick_dict_ppk.items():
        # read temp stream
        st_paths = sorted(glob.glob(os.path.join(temp_dir, '%s.*'%net_sta)))
        if len(st_paths)!=3: continue
        st = read_stream(st_paths)
        if len(st)!=3: continue
        # cut temp data
        temp_trig = trim_stream(st, tp-temp_win_trig[0], tp+temp_win_trig[1])
        temp_p = trim_stream(st, tp-temp_win_p[0], tp+temp_win_p[1])
        temp_s = trim_stream(st, ts-temp_win_s[0], ts+temp_win_s[1])
        temp = [st2np(sti) for sti in [temp_trig, temp_p, temp_s]]
        temp = [temp[i][:,0:temp_win_npts[i]] for i in range(3)]
        # calc norm
        norm_trig = np.array([sum(tr**2)**0.5 for tr in temp[0]])
        norm_p = np.array([sum(tr**2)**0.5 for tr in temp[1]])
        norm_s = np.array([sum(tr**2)**0.5 for tr in temp[2]])
        norm_temp = [norm_trig, norm_p, norm_s]
        # get time shift (dt)
        dt_list = [int(dt*samp_rate) for dt in [ot-tp+temp_win_trig[0], tp-ot, ts-ot]]
        pick_dict_data[net_sta] = [temp, norm_temp, dt_list]
    return temp_name, pick_dict_data

  def __len__(self):
    return len(self.temp_list)


""" Base functions
"""

# stream processing
def preprocess(stream):
    # time alignment
    start_time = max([trace.stats.starttime for trace in stream])
    end_time = min([trace.stats.endtime for trace in stream])
    if start_time>end_time: print('bad data!'); return []
    st = stream.slice(start_time, end_time)
    # resample data
    org_rate = int(st[0].stats.sampling_rate)
    rate = np.gcd(org_rate, samp_rate)
    if rate==1: print('warning: bad sampling rate!'); return []
    decim_factor = int(org_rate / rate)
    resamp_factor = int(samp_rate / rate)
    st = st.decimate(decim_factor)
    if resamp_factor!=1: st = st.interpolate(samp_rate)
    # filter
    st = st.detrend('demean').detrend('linear').taper(max_percentage=0.05, max_length=10.)
    flt_type, freq_rng = freq_band
    if flt_type=='highpass':
        return st.filter(flt_type, freq=freq_rng)
    if flt_type=='bandpass':
        return st.filter(flt_type, freqmin=freq_rng[0], freqmax=freq_rng[1])

def read_stream(stream_paths):
    stream  = read(stream_paths[0])
    stream += read(stream_paths[1])
    stream += read(stream_paths[2])
    return preprocess(stream)

def trim_stream(stream, start_time, end_time):
    return stream.copy().trim(start_time, end_time, pad=True, fill_value=0.)


# format transform
def cpu2cuda(data):
    return data.float().cuda(non_blocking=True)

def st2np(stream):
    return np.array([trace.data for trace in stream], dtype=np.float64)

