""" Configure file for hypoDD interface
"""
import os, sys
sys.path.append('/home/zhouyj/software/PAD')
# import PAD functions
import associators
import data_pipelines as dps

class Config(object):
  def __init__(self):

    # 1. mk_sta
    self.fsta_in = '/data2/Ridgecrest/station.csv'
    self.fsta_out = 'input/station.dat'

    # 2. mk_dt
    self.temp_pha = '../../run_pad/hypodd/output/rc_pad_reloc.pha'
    self.det_pha = '../output/rc/rc_msms.pha'
    self.date_rng = '20190704-20190711'
    self.num_proc = 7 # for parallel
    self.dep_corr = 5 # avoid air quake
    self.ot_dev = 2 # ot diff for det assoc
    self.cc_thres = 0.3 # select high cc event pairs
    self.dt_thres = [1., 1.5]
    self.max_nbr = 30
    self.fcsv = 'output/msms_reloc_fore.csv'
    sta_dict = dps.Data(None).get_rc_sta(self.fsta_in)
    self.associator = associators.TS_Assoc(sta_dict)
    if not os.path.exists('input'): os.makedirs('input')
    if not os.path.exists('output'): os.makedirs('output')

