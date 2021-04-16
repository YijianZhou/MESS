""" Configure file for hypoDD interface
"""
import os, sys
sys.path.append('/home/zhouyj/software/PAL')
import data_pipeline as dp

class Config(object):
  def __init__(self):

    # 1. mk_sta
    self.fsta_in = 'input/example_pal.sta'
    self.fsta_out = 'input/station.dat'

    # 2. mk_dt
    self.temp_pha = 'input/example_pal_hyp-ct_full.pha' # reloc template phase file
    self.det_pha = 'input/example_mess.pha' # mess output phase file
    self.time_range = '20190704-20190717'
    self.num_workers = 13
    self.evid_stride = 100000
    self.dep_corr = 5 # avoid air quake
    self.ot_dev = 2. # ot diff for det assoc
    self.cc_thres = 0.25 # select high cc event pairs
    self.dt_thres = [0.6,1.] # max dt_p & dt_s
    self.nbr_thres = [3,30] # min & max num of neighbor event
    self.min_sta = 4
    self.sta_dict = dp.get_sta_dict(self.fsta_in)

    # 3. reloc2csv
    self.ot_range = '20190704-20190717'
    self.lat_range = [35.4,36.1]
    self.lon_range = [-117.85,-117.25]
    self.xy_pad = [0.046,0.037] # degree
    self.num_grids = [25,25] # x,y (lon, lat)
    self.ctlg_code = 'example_mess_cc'
    self.keep_grids = False
