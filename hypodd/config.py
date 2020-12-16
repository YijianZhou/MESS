""" Configure file for hypoDD interface
"""
import os, sys
sys.path.append('/home/zhouyj/software/PAD')
# import PAD functions
import associators
import data_pipeline as dp

class Config(object):
  def __init__(self):

    # 1. mk_sta
    self.fsta_in = 'input/example_pad.sta'
    self.fsta_out = 'input/station.dat'

    # 2. mk_dt
    self.temp_pha = 'input/example_pad_hyp-ct_all.pha' # reloc template phase file
    self.det_pha = 'input/example_mess.pha' # mess output phase file
    self.time_range = '20190704-20190713'
    self.num_workers = 9
    self.evid_stride = 100000
    self.dep_corr = 5 # avoid air quake
    self.ot_dev = 2. # ot diff for det assoc
    self.cc_thres = 0.3 # select high cc event pairs
    self.dt_thres = [1.,1.8] # max dt_p & dt_s
    self.nbr_thres = [2,20] # min & max num of neighbor event
    sta_dict = dp.get_sta_dict(self.fsta_in)
    self.calc_mag = associators.TS_Assoc(sta_dict).calc_mag

    # 3. reloc2csv
    self.split_ranges = ['20190704-20190709','20190709-20190713']
    self.fctlg = ['output/example_mess_cc_0704-0708.ctlg',
        'output/example_mess_cc_0709-0712.ctlg'][0]
    self.write_temp = [False, True][1] # set one of the split as True

