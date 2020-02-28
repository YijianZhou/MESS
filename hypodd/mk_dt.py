from obspy import UTCDateTime
import numpy as np
import time
import argparse
import config

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--date_range', type=str,
                        default='20190704-20190705')
    parser.add_argument('--out_dt', type=str,
                        default='input/dt.cc')
    parser.add_argument('--out_event', type=str,
                        default='input/event.dat')
    parser.add_argument('--start_id', type=int, default=100000)
    args = parser.parse_args()


# i/o paths
cfg = config.Config()
dep_corr = cfg.dep_corr
out_dt = open(args.out_dt,'w')
out_event = open(args.out_event,'w')
# assoc params
ot_dev = cfg.ot_dev
cc_thres = cfg.cc_thres
dt_thres = cfg.dt_thres
nbr_thres = cfg.nbr_thres
start_date, end_date = [UTCDateTime(date) for date in args.date_range.split('-')]
associator = cfg.associator


# read temp pha
def mk_temp_dict():
    f=open(cfg.temp_pha); lines=f.readlines(); f.close()
    temp_dict = {}
    temp_id = 0
    for line in lines:
        codes = line.split(',')
        if len(codes)==6:
            temp_name = codes[0]
            ot = UTCDateTime(temp_name.split('_')[1])
            lat, lon, dep = [float(code) for code in codes[2:5]]
            temp_loc = [ot, lat, lon, dep]
            temp_dict[temp_name] = [temp_id, temp_loc, {}]
            temp_id += 1
        else:
            net_sta = '.'.join(codes[0:2])
            tp, ts = [UTCDateTime(code) for code in codes[2:4]]
            temp_dict[temp_name][-1][net_sta] = [tp, ts]
    return temp_dict


# read det pha
def mk_det_list():
    f=open(cfg.det_pha); lines=f.readlines(); f.close()
    dtype = [('temp','O'),('ot','O'),('loc','O'),('cc','O'),('picks','O')]
    det_list = []
    num=0
    for line in lines:
        codes = line.split(',')
        if len(codes[0])>16:
            temp_name = codes[0]
            if temp_name not in temp_dict: to_add = False; continue
            ot = UTCDateTime(codes[1])
            if ot < start_date: to_add = False; continue
            if ot > end_date: break
            to_add = True
            cc = float(codes[-1])
            det_list.append((temp_name, ot, temp_dict[temp_name][1][1:], cc, {}))
            num+=1
            if num%5000==0: print('processed %s dets'%num)
        else:
            if not to_add: continue
            net_sta = codes[0]
            tp, ts = [UTCDateTime(code) for code in codes[1:3]]
            s_amp, cc_p, cc_s = [float(code) for code in codes[3:6]]
            det_list[-1][-1][net_sta] = [tp, ts, s_amp, cc_p, cc_s]
    return np.array(det_list, dtype=dtype)


# write dt.cc
def write_dt(det, evid, fout):
    det_ot = det['ot']
    temp = temp_dict[det['temp']]
    temp_id, temp_ot = temp[0], temp[1][0]
    fout.write('# {:9} {:9} 0.0\n'.format(evid, temp_id))
    for net_sta, [tp, ts, _, cc_p, cc_s] in det['picks'].items():
        sta = net_sta.split('.')[1]
        temp_tp, temp_ts  = temp[-1][net_sta]
        temp_ttp, temp_tts = temp_tp-temp_ot, temp_ts-temp_ot
        det_ttp, det_tts = tp-det_ot, ts-det_ot
        dtp = det_ttp - temp_ttp
        dts = det_tts - temp_tts
        if abs(dtp)<dt_thres[0] and cc_p>cc_thres: 
            fout.write('{:7} {:8.5f} {:.4f} P\n'.format(sta, dtp, cc_p**0.5))
        if abs(dts)<dt_thres[1] and cc_s>cc_thres: 
            fout.write('{:7} {:8.5f} {:.4f} S\n'.format(sta, dts, cc_s**0.5))


# write event.dat
def write_event(event_loc, evid, fout):
    ot, lat, lon, dep, mag = event_loc
    dep += dep_corr
    date = '{:0>4}{:0>2}{:0>2}'.format(ot.year, ot.month, ot.day)
    time = '{:0>2}{:0>2}{:0>2}{:0>2.0f}'.format(ot.hour, ot.minute, ot.second, ot.microsecond/1e4)
    loc = '{:7.4f}   {:8.4f}   {:8.3f}  {:4.1f}'.format(lat, lon, dep, mag)
    err_rms = '   0.00    0.00   0.0'
    fout.write('{}  {}   {} {} {:>10}\n'.format(date, time, loc, err_rms, evid))


# calc mag with PAD assoc
def calc_mag(event):
    event_loc = {'evt_lat':event['loc'][0], 'evt_lon':event['loc'][1]}
    event_pick = [{'net_sta':net_sta, 's_amp':s_amp} \
        for net_sta, [_,_,s_amp,_,_] in event['picks'].items()]
    event_loc = associator.calc_mag(event_pick, event_loc)
    return event_loc['mag']


# prep input
print('prep input')
t=time.time()
temp_dict = mk_temp_dict()
det_list = mk_det_list()
det_list = det_list[[temp in temp_dict for temp in det_list['temp']]]
det_list = det_list[det_list['cc'] > cc_thres]

# assoc det
num_dets = len(det_list)
for i in range(num_dets):
    det_id = args.start_id + i
    if i%100==0: print('process {}th det | {:.0f} remains | {:.1f}s'\
        .format(i, num_dets, time.time()-t))
    det0 = det_list[0]

    # find nbr dets by ot cluster
    cond_ot = abs(det_list['ot'] - det0['ot']) < ot_dev
    reloc_dets = det_list[cond_ot]
    cc = reloc_dets['cc']
    cc_max = np.amax(cc)
    cc_min = np.sort(cc)[::-1][0:nbr_thres[1]][-1]
    deti = reloc_dets[np.argmax(cc)]
    reloc_dets = reloc_dets[cc >= cc_min]
    det_loc = [deti['ot']] + deti['loc'] + [calc_mag(deti)]

    # whether self-det
    temp_id, temp_loc = temp_dict[deti['temp']][0:2]
    is_self = abs(temp_loc[0]-det_loc[0])<ot_dev and cc_max>0.8
    if is_self:
        reloc_dets = reloc_dets[reloc_dets['cc'] < cc_max]
        det_id = temp_id

    # write dt.cc & event.dat
    if len(reloc_dets)>=nbr_thres[0] or is_self:
        for reloc_det in reloc_dets: write_dt(reloc_det, det_id, out_dt)
        write_event(det_loc, det_id, out_event)

    # next det
    det_list = np.delete(det_list, np.where(cond_ot), axis=0)
    if len(det_list)==0: break

out_dt.close()
out_event.close()
