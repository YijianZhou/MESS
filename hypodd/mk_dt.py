import sys
sys.path.append('/home/zhouyj/software/PAD')
from obspy import UTCDateTime
import numpy as np
import time
import argparse
import config
# import PAD functions
import associators
import data_pipeline as dp

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--temp_pha', type=str,
                        default='../../output/phase_temp_jz.dat')
    parser.add_argument('--det_pha', type=str,
                        default='../../output/phase_det_jz.dat')
    parser.add_argument('--time_range', type=str,
                        default='20170927,20170928')
    parser.add_argument('--out_dt', type=str,
                        default='input/dt.cc')
    parser.add_argument('--out_event', type=str,
                        default='input/event.dat')
    parser.add_argument('--start_id', type=int, default=0)
    args = parser.parse_args()


# i/o paths
cfg = config.Config()
dep_corr = cfg.dep_corr
out_dt = open(args.out_dt,'w')
out_event = open(args.out_event,'w')
# assoc params
ot_dev = cfg.ot_dev
cc_thres = cfg.cc_thres
start_date, end_date = [UTCDateTime(date) for date in args.time_range.split(',')]
# import associator for mag estimating
resp_dict = cfg.resp_dict
sta_dict = dp.get_jz_sta(cfg.fsta)
associator = associators.TS_Assoc(sta_dict, resp_dict)


# read temp pha
def mk_temp_dict():

  f=open(args.temp_pha); lines=f.readlines(); f.close()
  temp_dict = {}
  temp_id = 0
  for line in lines:
    codes = line.split(',')
    if len(codes)==5:
        temp_name = codes[0]
        ot = temp_name
        if ot[-6:]=='60.000': ot = ot[:-6]+'59.99'
        ot = UTCDateTime(ot)
        lat, lon, dep = [float(code) for code in codes[1:4]]
        temp_loc = [ot, lat, lon, dep]
        temp_dict[temp_name] = [temp_loc, temp_id, {}]
        temp_id += 1
    else:
        sta = codes[1]
        tp, ts = [UTCDateTime(code) for code in codes[2:4]]
        temp_dict[temp_name][-1][sta] = [tp, ts]
  return temp_dict


# read MSMS det pha
def mk_det_list():

  f=open(args.det_pha); lines=f.readlines(); f.close()
  dtype = [('temp','O'),('ot','O'),('lat','O'),('lon','O'),
           ('dep','O'),('cc','O'),('picks','O')]
  det_list = []
  num=0
  for line in lines:
    codes = line.split(',')
    if len(codes[0])>8:
        temp_name = codes[0]
        ot = UTCDateTime(codes[1])
        if ot < start_date: to_add = False; continue
        if ot > end_date: break
        to_add = True
        lat, lon, dep, cc = [float(code) for code in codes[2:6]]
        det_list.append((temp_name, ot, lat, lon, dep+dep_corr, cc, []))
        num+=1
        if num%5000==0: print('processed %s dets'%num)
    else:
        if not to_add: continue
        net, sta = codes[0].split('.')
        tp, ts = [UTCDateTime(code) for code in codes[1:3]]
        s_amp, cc_p, cc_s = [float(code) for code in codes[3:6]]
        pick_dict = {'net':net, 'sta':sta, 's_amp':s_amp, 
                     'tp':tp, 'ts':ts, 'cc_p':cc_p, 'cc_s':cc_s}
        det_list[-1][-1].append(pick_dict)
  return np.array(det_list, dtype=dtype)


# write dt.cc
def write_dt(det, evid, out):

    det_ot = det['ot']
    temp_name = det['temp']
    temp = temp_dict[temp_name]
    temp_ot, temp_id = temp[0][0], temp[1]
    out.write('# {:9} {:9} 0.0\n'.format(evid, temp_id))
    for pick_dict in det['picks']:
        sta = pick_dict['sta']
        temp_tp,  temp_ts  = temp[-1][sta]
        temp_ttp, temp_tts = temp_tp-temp_ot, temp_ts-temp_ot
        det_ttp = pick_dict['tp'] - det_ot
        det_tts = pick_dict['ts'] - det_ot
        cc_p, cc_s = pick_dict['cc_p'], pick_dict['cc_s']
        out.write('{:7} {:8.5f} {:.4f} P\n'.format(sta, det_ttp-temp_ttp, cc_p**0.5))
        out.write('{:7} {:8.5f} {:.4f} S\n'.format(sta, det_tts-temp_tts, cc_s**0.5))


# write event.dat
def write_event(event_loc, evid, out):

    ot, lat, lon, dep, mag = event_loc
    date = '{}{:0>2}{:0>2}'.format(ot.year, ot.month, ot.day)
    time = '{:0>2}{:0>2}{:0>2}{:0>2.0f}'.format(ot.hour, ot.minute, ot.second, ot.microsecond/1e4)
    err_rms = '   0.00    0.00   0.0'
    out.write('{}  {}   {:7.4f}   {:8.4f}   {:8.3f}  {:4.1f} {} {:>10}\n'\
    .format(date, time, lat, lon, dep, mag, err_rms, evid))


# calc mag with PAD assoc
def calc_mag(event):

    event_loc = {'evt_lon':event['lon'], 'evt_lat':event['lat']}
    event_pick = event[-1]
    event_loc = associator.calc_mag(event_pick, event_loc)
    return event_loc['mag']


# prep input
print('prep input')
t=time.time()
temp_dict = mk_temp_dict()
det_list = mk_det_list()

# assoc det
num_dets = len(det_list)
for i in range(num_dets):
    det_id = args.start_id + i
    if i%100==0: print('process {}th det | {:.0f} remains | {:.1f}s'\
               .format(i, num_dets, time.time()-t))
    det0 = det_list[0]

    # find nbr dets by ot cluster
    cond_ot = abs(det_list['ot'] - det0['ot']) < ot_dev
    cc = det_list[cond_ot]['cc']
    cc_max = np.amax(cc)
    reloc_dets = det_list[cond_ot]
    deti = reloc_dets[np.argmax(cc)]
    det_loc = [deti['ot'], deti['lat'], deti['lon'], deti['dep'], calc_mag(deti)]

    # whether self-det
    is_self = False
    temp_loc, temp_id = temp_dict[deti['temp']][0:2]
    if abs(temp_loc[0]-det_loc[0])<ot_dev and cc_max>0.9: is_self = True
    if is_self: 
        reloc_dets = np.delete(reloc_dets, np.argmax(cc), axis=0)
        det_id = temp_id

    # write dt.cc & event.dat
    for reloc_det in reloc_dets:
      if reloc_det['cc'] > cc_thres:
        write_dt(reloc_det, det_id, out_dt)
    if cc_max>cc_thres: write_event(det_loc, det_id, out_event)

    # next det
    det_list = np.delete(det_list, np.where(cond_ot), axis=0)
    if len(det_list)==0: break

out_dt.close()
out_event.close()
