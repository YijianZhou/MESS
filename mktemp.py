import os, sys, glob, shutil
sys.path.append('/home/zhouyj/software/mycode')
sys.path.append('/home/zhouyj/Xiaojiang/PpkAssocLoc')
import numpy as np
from obspy.core import *
import sac
import pickers
import data_pipeline as dp

def calc_tt(dist):
    vp = 5.9
    vs = 3.4
    ttp = dist/vp
    tts = dist/vs
    return ttp, tts

def calc_ot(tp, ts):
    vp = 5.9
    vs = 3.4
    r = (ts-tp) /(1/vs - 1/vp)
    ttp = r/vp
    return tp-ttp


# i/o files
net = 'XLS'
catalog = './input/catalog_%s.dat'%net
phase   = './input/phase_%s.dat'%net
f=open(catalog); events=f.readlines(); f.close()
f=open(phase);   phases=f.readlines(); f.close()
out_dir = './output/%s'%net
temp_dir = os.path.join(out_dir, 'Templates')
if not os.path.exists(temp_dir):
    os.makedirs(temp_dir)
ctlg_path = os.path.join(out_dir, 'ctlg_%s_temp.dat'%net)
pha_path  = os.path.join(out_dir, 'pha_%s_temp.dat'%net)
out_ctlg = open(ctlg_path, 'w')
out_pha  = open(pha_path, 'w')
# make sta dir
sta_file = '/data3/XJ_SAC/header/station_%s.dat'%net
sta_dict = dp.get_sta_dict(sta_file)

# params
max_dist = 100 # km
win_len = 45 # sec
min_sta = 4 # min number of sta
p_dt = 15 # sec before tp (cut win)
pick_win = [10*100, 1*100] # [lwin, swin], in points
p_win = [150, 150] # in points
s_win = [200, 200] # in points
p_thres = 5
s_thres = 3
pick_thres = 0.97
s_stride = 2
data_dir = '/data3/XJ_SAC/*/*'
picker = pickers.Trad_PS()

# for all events
event_id=0
event_picks=[]
for phase in phases:

  # head line
  if len(phase.split(','))==2:
    num_pick = int(phase.split(',')[1]); i=0
    event = events[event_id]
    event_id+=1
    print('-'*40)
    print('checking event {}'.format(event[:-1]))
    ot, lat, lon, dep = event.split(',')
    event_dir = os.path.join(temp_dir, ot)
    if not os.path.exists(event_dir): os.mkdir(event_dir)
    event_picks.append([event, event_dir, []])#TODO
    # convert format
    ot =  UTCDateTime(ot)
    lat = float(lat)
    lon = float(lon)
    dep = float(dep)
    data_dict = dp.get_xj(data_dir, ot)

  # org picks
  else:
    # select all picks
    sta = phase.split(',')[1]
    sta_line = sta_dict[sta_dict['station']==sta]
    sta_lat = float(sta_line['latitude'])
    sta_lon = float(sta_line['longitude'])
    sta_ele = float(sta_line['elevation'])
    dist0 = 111*np.sqrt((sta_lat-lat)**2 + (sta_lon-lon)**2)
    if dist0>max_dist: continue

    # theoretical arrival time
    dist1 = np.sqrt(dist0**2 + (dep+2)**2)
    ttp, tts = calc_tt(dist1) # travel time
    tp0 = ot + ttp # arrival time
    ts0 = ot + tts
    b = tp0 -UTCDateTime(ot.date) -p_dt
    # cut event
    out_paths=[]
    for data in data_dict[sta]:
        data = os.path.split(data)[-1]
        net, chn = data.split('.')[0], data.split('.')[-2]
        out_paths.append(os.path.join(event_dir,'%s.%s.%s'%(net,sta,chn)))
    sac.cut(data_dict[sta][0], b, b+win_len, out_paths[0])
    sac.cut(data_dict[sta][1], b, b+win_len, out_paths[1])
    sac.cut(data_dict[sta][2], b, b+win_len, out_paths[2])
    # read event data
    st = read(out_paths[0])
    st+= read(out_paths[1])
    st+= read(out_paths[2])
    st.detrend('constant').filter('highpass', freq=1.)

    # calc P snr
    p_data = st[2].data[p_dt*100 -p_win[0] -pick_win[0]
                       :p_dt*100 +p_win[1] +pick_win[1]]
    p_cf  = picker.calc_cf(p_data, pick_win)
    p_snr = np.amax(p_cf)
    p_idx = -pick_win[0] -p_win[0] +\
            np.where(p_cf >= pick_thres*np.amax(p_cf))[0][0]
    tp = tp0 + p_idx/100

    # calc S snr
    p_idx = p_idx + p_dt*100 # idx to time_win start_time
    picker.s_win = [int(100*(tp-ts0)+s_win[0]), int(100*(ts0-tp)+s_win[1])]
    st_data = [st[0].data, st[1].data, st[2].data]
    s_data0, s_data1 = picker.calc_filter(st_data, p_idx, s_stride)
    s_cf0 = picker.calc_cf(s_data0, pick_win, decim=s_stride)
    s_cf1 = picker.calc_cf(s_data1, pick_win, decim=s_stride)
    s_snr0 = np.amax(s_cf0)
    s_snr1 = np.amax(s_cf1)
    s_idx = -pick_win[0] -s_win[0] + s_stride*int(0.5*(\
            np.where(s_cf0 >= pick_thres *np.amax(s_cf0))[0][0] +\
            np.where(s_cf1 >= pick_thres *np.amax(s_cf1))[0][0]))
    ts = ts0 + s_idx/100
    if ts-tp<1: ts=ts0

    # output phase line
    print('{}: P snr {:.1f} dt {:.2f}; S snr {:.1f} {:.1f} dt {:.2f}'.\
         format(sta, p_snr, tp-tp0, s_snr0, s_snr1, ts-ts0))
    t0 = tp-tp0 +p_dt
    t1 = ts-tp0 +p_dt
    sac.ch_event(out_paths[0], lon, lat, dep, 1, t0, t1)
    sac.ch_event(out_paths[1], lon, lat, dep, 1, t0, t1)
    sac.ch_event(out_paths[2], lon, lat, dep, 1, t0, t1)
    oti = calc_ot(tp, ts)
    if  p_snr>p_thres\
    and s_snr0>s_thres\
    and s_snr1>s_thres:
        event_picks[-1][-1].append('{},{},{},{},{},{},{:.1f},{:.1f},{:.1f}\n'\
            .format(sta, tp, ts, tp0, ts0, oti, p_snr, s_snr0, s_snr1))
    else:
        for fname in out_paths: os.unlink(fname)


# select all event picks
for event_pick in event_picks:
    event = event_pick[0]
    event_dir = event_pick[1]
    picks = event_pick[2]
    if len(picks)<min_sta:
        shutil.rmtree(event_dir)
        continue
    out_ctlg.write(event)
    out_pha.write(event)
    for pick in picks:
        out_pha.write(pick)

