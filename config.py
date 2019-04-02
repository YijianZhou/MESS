class Config(object):
  def __init__(self):
    self.dt_trig = [1., 9.]
    self.dt_p = [1.,2.]
    self.dt_s = [0, 3]
    self.trig_thres = 0.35
    self.mask_len = 1.
    # data process
    self.decim_rate = 2
    self.freq_band = [1., 40.]

