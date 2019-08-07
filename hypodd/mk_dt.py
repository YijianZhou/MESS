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
    parser.add_argument('--out_ctlg', type=str,
                        default='event.sel')
    parser.add_argument('--out_dt', type=str,
                        default='dt.cc')
    parser.add_argument('--write_temp', type=int, default=1)
    parser.add_argument('--start_id', type=int, default=0)
    args = parser.parse_args()


# i/o paths
out_ctlg = open(args.out_ctlg,'w')
out_dt = open(args.out_dt,'w')

# assoc params
ot_dev = 5 # sec
loc_dev = 15 # km
dep_dev = 20 # km
start_date = UTCDateTime(args.time_range.split(',')[0])
end_date   = UTCDateTime(args.time_range.split(',')[1])

# import detector for mag estimating
resp_dict = {'CI':1e2, 'NN':1e2}    
sta_dict = dp.get_ci_sta(args.sta_file)
detector = detectors.TS_Det(sta_dict, resp_dict)


""" Prep MFT output files
"""
def mk_temp_dict():
  f=open(args.temp_pha); lines=f.readlines(); f.close()
  dtype=[('evid','O'),('ot','O'),('lat','O'),('lon','O'),('dep','O')]
  temp_dict = {}
  evid = 1
  for line in lines:
    info = line.split(',')
    if len(info)==5:
        temp_name = info[0]
        ot = UTCDateTime(info[0])
        lat = float(info[1])
        lon = float(info[2])
        dep = float(info[3])
        temp_loc = np.array([(evid, ot, lat, lon, dep)], dtype=dtype)
        temp_dict[temp_name] = [temp_loc, {}]
        evid += 1
        if evid%5000==0: print('processed %s temps'%evid)
    else:
        sta = info[0]
        tp = UTCDateTime(info[1])
        ts = UTCDateTime(info[2])
        temp_dict[temp_name][1][sta] = [tp, ts]
  return temp_dict


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
def write_ctlg(event, evid, out, get_mag=False):
    ot   = event['ot']
    lat  = event['lat']
    lon  = event['lon']
    dep  = event['dep']
    yr  = ot.year
    mag = -1.0
    if get_mag: mag = calc_mag(event)
    mon = str(ot.month).zfill(2)
    day = str(ot.day).zfill(2)
    date = '{}{:0>2}{:0>2}'.format(yr, mon, day)
    hr  = ot.hour
    mn  = ot.minute
    sec = ot.second + ot.microsecond/1e6
    time = '{:0>2}{:0>2}{:0>2}{:0>2}'.format(hr, mn, int(sec), int(100*sec%100))
    out.write('{}  {}  {:8.4f}  {:9.4f}  {:9.3f}  {:4.1f}    0.00    0.00   0.00  {:9}\n'\
      .format(date, time, lat, lon, dep, mag, evid))


def write_dt(det, temp_dict, evid, out):
    det_ot = det['ot']
    temp = temp_dict[det['temp']]
    temp_ot = temp[0][0]['ot']
    out.write('# {:9} {:9} 0.0\n'.format(evid, temp[0][0]['evid']))
    for pick_dict in det['picks']:
        sta = pick_dict['sta']
        [temp_tp, temp_ts] = temp[1][sta]
        temp_ttp, temp_tts = temp_tp-temp_ot, temp_ts-temp_ot
        det_ttp = pick_dict['tp'] - det_ot
        det_tts = pick_dict['ts'] - det_ot
        cc_p, cc_s = pick_dict['cc_p'], pick_dict['cc_s']
        out.write('{:7} {:8.5f} {:.4f} P\n'.format(sta, det_ttp-temp_ttp, cc_p**0.5))
        out.write('{:7} {:8.5f} {:.4f} S\n'.format(sta, det_tts-temp_tts, cc_s**0.5))


def calc_mag(event):
    event_loc = {'evt_lon':event['lon'], 'evt_lat':event['lat']}
    event_pick = event[-1]
    event_loc = detector.calc_mag(event_pick, event_loc)
    return event_loc['mag']


# prep input
print('prep input')
t=time.time()
temp_dict = mk_temp_dict()
det_list = mk_det_list()

# write temp events
if args.write_temp==1:
  print('write temp events: event.sel')
  for i,temp_name in enumerate(temp_dict):
    write_ctlg(temp_dict[temp_name][0][0], i+1, out_ctlg)

# assoc det
evid = args.start_id
for i,_ in enumerate(det_list):
    if i%100==0: print('process {}th det | {:.0f} remains | {:.1f}s'\
               .format(i, len(det_list), time.time()-t))
    det0 = det_list[0]
    # find nbr dets
    cond_ot = abs(det_list['ot'] - det0['ot']) < ot_dev
    cc = det_list[cond_ot]['cc']
    detsi = det_list[cond_ot]
    deti  = detsi[np.argmax(cc)]
    write_ctlg(deti, evid, out_ctlg, get_mag=True)
    # find dets for reloc
    reloc_dets = det_list[cond_ot]
    cond_loc = 111*((reloc_dets['lat'] - deti['lat'])**2 \
                  + (reloc_dets['lon'] - deti['lon'])**2)**0.5 < loc_dev
    cond_dep = abs(reloc_dets['dep'] - deti['dep']) < dep_dev
    reloc_dets = reloc_dets[cond_loc * cond_dep]
    for reloc_det in reloc_dets: write_dt(reloc_det, temp_dict, evid, out_dt)
    # next det
    evid += 1
    det_list = np.delete(det_list, np.where(cond_ot), axis=0)
    if len(det_list)==0: break

out_ctlg.close()
out_dt.close()
