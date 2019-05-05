import os, sys, glob, shutil
sys.path.append('/home/zhouyj/software/mycode')
sys.path.append('/home/zhouyj/Xiaojiang/PpkAssocLoc')
import argparse
import numpy as np
from obspy.core import read, UTCDateTime
import sac
import pickers
import data_pipeline as dp
import config

def calc_tt(dist):
    vp = 5.9
    vs = 3.4
    ttp = dist/vp
    tts = dist/vs
    return ttp, tts


def main(args):

    # i/o files
    data_dir = args.data_dir
    org_ctlg = args.org_ctlg
    f=open(org_ctlg); events=f.readlines(); f.close()
    out_root = args.out_root
    temp_root = os.path.join(out_root, args.temp_root)
    if not os.path.exists(temp_root): os.makedirs(temp_root)
    ctlg_path = os.path.join(out_root, args.out_ctlg)
    pha_path  = os.path.join(out_root, args.out_pha)
    out_ctlg = open(ctlg_path,'w')
    out_pha  = open(pha_path,'w')
    # make sta dir
    sta_file = args.sta_file
    sta_dict = dp.get_sta_dict(sta_file)

    # params
    cfg = config.Config()
    max_dist = cfg.max_dist # km
    win_len = cfg.win_len # sec
    min_sta = cfg.min_sta # min number of sta
    t_blank = cfg.t_blank # sec before tp (cut win)
    dt_p = cfg.dt_s
    dt_s = cfg.dt_s
    p_thres = cfg.p_thres
    s_thres = cfg.s_thres
    fd_thres = cfg.fd_thres
    picker = pickers.Trad_PS(trig_thres = p_thres, 
                             fd_thres = fd_thres)

    # for all events
    event_picks=[]
    for i, event in enumerate(events):

      # event info
      ot, lat, lon, dep, mag = event.split(',')
      event_dir = os.path.join(temp_root, ot)
      if not os.path.exists(event_dir): os.mkdir(event_dir)
      ot = UTCDateTime(ot)
      lat = float(lat)
      lon = float(lon)
      dep = float(dep)
      data_dict = dp.get_xj(data_dir, ot)
      event_picks.append([event, event_dir, []])
      print('-'*40)
      print('checking {}th event'.format(i))

      for sta in data_dict:
        # select all stations
        sta_line = sta_dict[sta_dict['station']==sta]
        sta_lat = float(sta_line['latitude'])
        sta_lon = float(sta_line['longitude'])
        sta_ele = float(sta_line['elevation'])
        dist0 = 111*np.sqrt((sta_lat-lat)**2 + (sta_lon-lon)**2)
        if dist0>max_dist: continue

        # theoretical arrival time
        dist1 = np.sqrt(dist0**2 + (dep+2.5)**2)
        ttp, tts = calc_tt(dist1) # travel time
        tp0 = ot + ttp # arrival time
        ts0 = ot + tts
        b = tp0 - UTCDateTime(ot.date) - t_blank
        # cut event
        out_paths = []
        for data in data_dict[sta]:
            data = os.path.split(data)[-1]
            net, chn = data.split('.')[0], data.split('.')[-2]
            out_paths.append(os.path.join(event_dir,'%s.%s.%s'%(net,sta,chn)))
        if len(data_dict[sta])!=3: continue
        sac.cut(data_dict[sta][0], b, b+win_len, out_paths[0])
        sac.cut(data_dict[sta][1], b, b+win_len, out_paths[1])
        sac.cut(data_dict[sta][2], b, b+win_len, out_paths[2])
        # read event data
        st = read(out_paths[0])
        st+= read(out_paths[1])
        st+= read(out_paths[2])

        # run picker
        picks = picker.pick(st)
        if len(picks)==0: 
            for fname in out_paths: os.unlink(fname)
            continue
        tp = picks['p_arr'][0]
        ts = picks['s_arr'][-1]
        p_snr = picks['p_snr'][0]
        s_snr = picks['s_snr'][-1]

        # write head
        t0 = tp - tp0 + t_blank
        t1 = ts - tp0 + t_blank
        sac.ch_event(out_paths[0], lon, lat, dep, 1, [t0,t1])
        sac.ch_event(out_paths[1], lon, lat, dep, 1, [t0,t1])
        sac.ch_event(out_paths[2], lon, lat, dep, 1, [t0,t1])

        # select
        if s_snr>s_thres\
        and abs(tp-tp0)<dt_p\
        and abs(ts-ts0)<dt_s:
            event_picks[-1][-1].append('{},{},{},{:.1f},{:.1f}\n'\
                .format(sta, tp, ts, p_snr, s_snr))
        else:
            for fname in out_paths: os.unlink(fname)

    # select all event picks
    for event_pick in event_picks:
        print('select {}'.format(event_pick[0][:-1]))
        event = event_pick[0]
        event_dir = event_pick[1]
        picks = event_pick[2]
        if len(picks)<min_sta:
          if os.path.exists(event_dir):
            shutil.rmtree(event_dir)
          continue
        out_ctlg.write(event)
        out_pha.write(event)
        for pick in picks: out_pha.write(pick)

    out_ctlg.close()
    out_pha.close()


if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--data_dir', type=str,
                        default='/data3/XJ_SAC/[Y-Z]*/*')
    parser.add_argument('--org_ctlg', type=str,
                        default='/home/zhouyj/Xiaojiang/PpkAssocLoc/hypoinverse/xj_zsy.csv')
    parser.add_argument('--sta_file', type=str,
                        default='/data3/XJ_SAC/header/station_ZSY.dat')
    parser.add_argument('--out_root', type=str,
                        default='./output')
    parser.add_argument('--temp_root', type=str,
                        default='Templates')
    parser.add_argument('--out_ctlg', type=str,
                        default='zsy_temp.ctlg')
    parser.add_argument('--out_pha', type=str,
                        default='zsy_temp.pha')
    args = parser.parse_args()
    main(args)

