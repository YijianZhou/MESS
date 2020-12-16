import sys
sys.path.append('/home/zhouyj/software/PAD')
import data_pipeline as dp
import pickers

class Config(object):
  def __init__(self):

    # Template cut
    self.win_len = [15,25]  # cut window length
    self.min_sta = 4        # min sta num for a template events
    self.max_sta = 15
    self.get_data_dict = dp.get_data_dict
    self.get_sta_dict = dp.get_sta_dict

    # MESS params
    self.temp_win_det = [1.,11.]   # win for detection, pre & post P
    self.temp_win_p = [0.5,1.5]    # win for p pick, pre & post P
    self.temp_win_s = [0.5,2.5]    # win for s pick, pre & post S
    self.trig_thres = 0.25         # cc thres for det & peak expansion
    self.expand_len = 2.5          # win len for cc peak expansion
    self.det_gap = 5.              # gap sec for detection
    self.pick_win_p = [1.5, 1.5]   # win for P pick
    self.pick_win_s = [2.5, 2.5]   # win for S pick
    self.chn_p = [2]
    self.chn_s = [0,1]
    self.amp_win = [1, 5]          # win for get_amp, sec pre&post P

    # data process
    self.samp_rate = 50
    self.freq_band = ['bandpass', [1., 40.]]
    self.get_amp = pickers.STA_LTA_PCA().get_amp
    self.num_workers = 10

