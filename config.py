class Config(object):
  def __init__(self):

    # Template selection params
    self.win_len  = 45       # cut window length
    self.t_blank  = 15       # time before P in the cut window
    self.min_sta  = 4        # min sta num for a template events

    # MFT params
    self.win_trig = [0.5, 5.]
    self.win_p = [0.5,1.]
    self.win_s = [0, 1.5]
    self.trig_thres = 0.25
    self.mask_len = 0.5

    # data process
    self.decim_rate = 2
    self.freq_band = ['bandpass', 1., 40.]

