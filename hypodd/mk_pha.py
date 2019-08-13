import sys
sys.path.append('/home/zhouyj/Documents/PAD')
from obspy import UTCDateTime
import numpy as np
import time
import argparse
# import PAD
import detectors
import data_pipeline as dp

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--temp_pha', type=str,
                        default='/home/zhouyj/Desktop/California/MFT/output/ci/ci_temp.pha')
    parser.add_argument('--det_pha', type=str,
                        default='/home/zhouyj/Desktop/California/MFT/output/ci2/phase_aug_rc.dat')
    parser.add_argument('--sta_file', type=str,
                        default='/home/zhouyj/Desktop/California/preprocess/station.dat')
    parser.add_argument('--time_range', type=str,
                        default='20190704,20190713')
    parser.add_argument('--out_pha', type=str,
                        default='rc.pha')
    parser.add_argument('--out_temp', type=str,
                        default='input2/temp.pha')
    parser.add_argument('--write_temp', type=int, default=1)
    parser.add_argument('--start_id', type=int, default=0)
    args = parser.parse_args()


# i/o paths
out_pha = open(args.out_pha,'w')

# assoc params
ot_dev = 5 # sec
cc_thres = 0.23
start_date = UTCDateTime(args.time_range.split(',')[0])
end_date   = UTCDateTime(args.time_range.split(',')[1])

# import detector for mag estimating
resp_dict = {'CI':1e2, 'NN':1e2}    
sta_dict = dp.get_ci_sta(args.sta_file)
detector = detectors.TS_Det(sta_dict, resp_dict)


""" Prep MFT output files
"""
def mk_temp_list():
  f=open(args.temp_pha); lines=f.readlines(); f.close()
  dtype=[('ot','O'),('lat','O'),('lon','O'),('dep','O'),('picks','O')]
  temp_list = []
  num=0
  for line in lines:
    info = line.split(',')
    if len(info)==5:
        temp_name = info[0]
        ot = UTCDateTime(info[0])
        lat = float(info[1])
        lon = float(info[2])
        dep = float(info[3])
        temp_list.append((ot, lat, lon, dep, []))
        num+=1
        if num%5000==0: print('processed %s temps'%num)
    else:
        sta = info[0]
        tp = UTCDateTime(info[1])
        ts = UTCDateTime(info[2])
        pick_dict = {'sta':sta, 'tp':tp, 'ts':ts}
        temp_list[-1][-1].append(pick_dict)
  return np.array(temp_list, dtype=dtype)


def mk_det_list():
  f=open(args.det_pha); lines=f.readlines(); f.close()
  dtype = [('temp','O'),('ot','O'),('lat','O'),('lon','O'),
           ('dep','O'),('cc','O'),('picks','O')]
  det_list = []
  num=0
  for line in lines:
    info = line.split(',')
    if len(info[0])>8:
        temp_name = info[0]
        ot = UTCDateTime(info[1])
        if ot < start_date: to_add = False; continue
        if ot > end_date: break
        to_add = True
        lat = float(info[2])
        lon = float(info[3])
        dep = float(info[4])
        cc  = float(info[5])
        det_list.append((temp_name, ot, lat, lon, dep, cc, []))
        num+=1
        if num%5000==0: print('processed %s dets'%num)
    else:
        if not to_add: continue
        sta = info[0] #TODO
        net = 'CI'
        if sta in ['GWY','QSM']: net='NN'
        tp = UTCDateTime(info[1])
        ts = UTCDateTime(info[2])
        s_amp = float(info[3])
        cc_p = float(info[4])
        cc_s = float(info[5])
        pick_dict = {'net':net, 'sta':sta, 's_amp':s_amp, 
                     'tp':tp, 'ts':ts, 'cc_p':cc_p, 'cc_s':cc_s}
        det_list[-1][-1].append(pick_dict)
  return np.array(det_list, dtype=dtype)


""" Write hypoDD input files
"""
def write_pha(event, evid, out, get_mag=False):

    # 1. write event info
    ot   = event['ot']
    lat  = event['lat']
    lon  = event['lon']
    dep  = event['dep']
    yr  = ot.year
    mag = -1.0
    if get_mag: mag = calc_mag(event)
    mon = ot.month
    day = ot.day
    date = '{} {:2} {:2}'.format(yr, mon, day)
    hr  = ot.hour
    mn  = ot.minute
    sec = ot.second + ot.microsecond/1e6
    time = '{:2} {:2} {:5.2f}'.format(hr, mn, sec)
    out.write('# {} {}  {:7.4f} {:9.4f}  {:6.2f} {:4.2f}  0.00  0.00  0.00  {:>9}\n'\
      .format(date, time, lat, lon, dep, mag, evid))

    # 2. write pha info
    for pick_dict in event['picks']:
        sta = pick_dict['sta']
        ttp = pick_dict['tp'] - ot
        tts = pick_dict['ts'] - ot
        out.write('{:<5}{}{:6.3f}  {:6.3f}   P\n'.format(sta, ' '*6, ttp, 1.))
        out.write('{:<5}{}{:6.3f}  {:6.3f}   S\n'.format(sta, ' '*6, tts, 1.))


def calc_mag(event):
    event_loc = {'evt_lon':event['lon'], 'evt_lat':event['lat']}
    event_pick = event[-1]
    event_loc = detector.calc_mag(event_pick, event_loc)
    return event_loc['mag']


# prep input
print('prep input')
t=time.time()
det_list = mk_det_list()

# write temp events
if args.write_temp==1:
  print('write temp events: event.sel')
  out_temp = open(args.out_temp, 'w')
  temp_list = mk_temp_list()
  for i,temp in enumerate(temp_list):
    write_pha(temp, i+1, out_temp)
  out_temp.close()

# assoc det
evid = args.start_id
for i,_ in enumerate(det_list):
    if i%100==0: print('process {}th det | {:.0f} remains | {:.1f}s'\
               .format(i, len(det_list), time.time()-t))
    det0 = det_list[0]
    # find nbr dets
    cond_ot = abs(det_list['ot'] - det0['ot']) < ot_dev
    cc = det_list[cond_ot]['cc']
    cc_max = np.amax(cc)
    detsi = det_list[cond_ot]
    deti  = detsi[np.argmax(cc)]
    if cc_max>cc_thres:
        write_pha(deti, evid, out_pha, get_mag=True)
    # next det
    evid += 1
    det_list = np.delete(det_list, np.where(cond_ot), axis=0)
    if len(det_list)==0: break

out_pha.close()

