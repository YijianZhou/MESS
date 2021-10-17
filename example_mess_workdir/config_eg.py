import sys
sys.path.append('/home/zhouyj/software/PAL')
import data_pipeline as dp

class Config(object):
  def __init__(self):

    # Template cut
    self.win_len = [10,20]  # cut window length
    self.min_sta = 4        # min sta num for a template events
    self.max_sta = 15
    self.get_data_dict = dp.get_data_dict
    self.get_sta_dict = dp.get_sta_dict

    # MFT params
    self.temp_win_det = [1.,9.]    # win for detection, pre & post P
    self.temp_win_p = [0.5,1.5]    # win for p pick, pre & post P
    self.temp_win_s = [0.5,2.5]    # win for s pick, pre & post S
    self.trig_thres = 0.3          # cc thres for det & peak expansion
    self.expand_len = 1.           # win len for cc peak expansion
    self.det_gap = 5.              # gap sec for detection
    self.pick_win_p = [1.0,1.0]    # win for P pick
    self.pick_win_s = [1.6,1.6]    # win for S pick
    self.chn_p = [2]
    self.chn_s = [0,1]
    self.amp_win = [1,4]

    # data process
    self.to_prep = False  # if to prep template in run_mess stage
    self.samp_rate = 50
    self.freq_band = [2.,40.]
    self.num_workers = 10

