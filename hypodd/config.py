class Config(object):
  def __init__(self):

    self.dep_corr = 0 # avoid air quake
    self.ot_dev = 3 # ot diff for det assoc
    self.cc_thres = 0.35 # select high cc event pairs
    self.resp_dict = {'JZ':1e5}
    self.fsta = 'input/station.dat'
    self.fcsv = 'output/jz_msms_reloc.csv'

