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
    self.fsta_in = '/data2/ZSY_SAC/header/station_ZSY.dat'
    self.fsta_out = 'input/station.dat'

    # 2. mk_dt
    self.temp_pha = '../../run_pad/hypodd/output/zsy_reloc.pha'
    self.det_pha = '../output/zsy/catalog_20180325-20190201.dat'
    self.dep_corr = 5 # avoid air quake
    self.ot_dev = 3 # ot diff for det assoc
    self.cc_thres = 0.3 # select high cc event pairs
    self.fcsv = 'output/zsy_msms_reloc.csv'
    dp = dps.Data(None)
    sta_dict = dp.get_sta_dict(self.fsta_in)
    self.associator = associators.TS_Assoc(sta_dict)
    if not os.path.exists('input'): os.makedirs('input')
    if not os.path.exists('output'): os.makedirs('output')

