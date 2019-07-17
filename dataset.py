import os, glob, time
from multiprocessing import Pool
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


def st2np(stream):
    return np.array([stream[0].data, stream[1].data, stream[2].data])


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

def read_temp(temp_dir, temp_dict):
    temp_name = os.path.split(temp_dir)[-1]
    out_dict = {}
    # read data
    temp_paths = glob.glob(os.path.join(temp_dir,'*.*Z'))
    for temp_path in temp_paths:
        net, sta, _ = os.path.split(temp_path)[-1].split('.')
        [tp, ts, dt_ot] = temp_dict[temp_name][1][sta]
        stream_paths = sorted(glob.glob(os.path.join(temp_dir, '%s.%s.*'%(net,sta))))
        stream = read(stream_paths[0])
        stream+= read(stream_paths[1])
        stream+= read(stream_paths[2])
        stream = preprocess(stream)
        temp_trig = stream.slice(tp - win_trig[0], tp + win_trig[1])
        temp_p    = stream.slice(tp - win_p[0], tp + win_p[1])
        temp_s    = stream.slice(ts - win_s[0], ts + win_s[1])
        if len(temp_trig)!=3 or len(temp_s)!=3: continue
        temp = [temp_trig, temp_p, temp_s]
        temp = [st2np(tempi) for tempi in temp]
        norm_trig = [np.sqrt(np.sum(tri**2)) for tri in temp[0]]
        norm_p    = [np.sqrt(np.sum(tri**2)) for tri in temp[1]]
        norm_s    = [np.sqrt(np.sum(tri**2)) for tri in temp[2]]
        norm_temp = [norm_trig, norm_p, norm_s]
        # output (data, label)
        out_dict[sta] = [temp, norm_temp]
    return [temp_name, out_dict]


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
    pool = Pool(processes=10)
    temp_dirs = glob.glob(os.path.join(temp_root,'*'))
    num_temp = len(temp_dirs)
    print('read template data')
    t=time.time()
    out = []
    for temp_dir in temp_dirs:
        outi = pool.apply_async(read_temp, args=(temp_dir, temp_dict))
        out.append(outi)
    pool.close()
    pool.join()
    print('Read {} templates | time comsumption: {:.1f}s'.format(num_temp, time.time()-t))

    # unpack results
    for outi in out:
        [temp_name, out_dict] = outi.get()
        [ot, lat, lon, dep] = temp_dict[temp_name][0]
        for sta in out_dict:
            [tp, ts, dt_ot] = temp_dict[temp_name][1][sta]
            out_dict[sta] += [tp, ts, ot, dt_ot]
        temp_dict[temp_name][1] = out_dict

    return temp_dict


""" Read Continuous Data
  Inputs
    data_dict: data_dict[sta] = stream_paths
"""

def get_sta_data(data_paths, date):
    stream = read(data_paths[0])
    stream+= read(data_paths[1])
    stream+= read(data_paths[2])
    stream = preprocess(stream)
    st_np = st2np(stream)
    data_cum = [np.cumsum(datai**2) for datai in st_np]
    norm_data = [np.sqrt(cumi[npts_trig:] - cumi[:-npts_trig]) for cumi in data_cum]
    dt_st = stream[0].stats.starttime - date
    return [stream, norm_data, date, dt_st]


# read continuous raw data and preprocess
def read_data(data_dict, date):

    print('read continuous data')
    t=time.time()
    num_sta = len(data_dict)
    pool = Pool(processes = num_sta)
    out = []
    for sta in data_dict:
        outi = pool.apply_async(get_sta_data, args=(data_dict[sta],date))
        out.append([sta, outi])
    pool.close()
    pool.join()
    for [sta, outi] in out: data_dict[sta] = outi.get()
    print('Read {} stations | time comsumption: {:.1f}s'.format(num_sta, time.time()-t))
    todel = [sta for sta in data_dict if len(data_dict[sta][0])!=3]
    for sta in todel: data_dict.pop(sta)
    return data_dict

