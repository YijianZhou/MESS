class Config(object):
  def __init__(self):
    # Template selection params
    self.max_dist = 100      # max epicentral distance for sta obs
    self.win_len  = 45       # cut window length
    self.t_blank  = 15       # time before P in the cut window
    self.min_sta  = 4        # min sta num for a template events
    self.p_thres  = 20       # min P STA/LTA (=trig_thres in PAL)
    self.s_thres  = 10       # min S STA/LTA
    self.fd_thres = 3.       # min value of dominant frequency
    self.dt_p     = 1.5      # P arrival misfit
    self.dt_s     = 3.       # S arrival misfit
    # MFT params
    self.dt_trig = [1., 9.]
    self.dt_p = [1.,2.]
    self.dt_s = [0, 3]
    self.trig_thres = 0.3
    self.mask_len = 1.

    # data process
    self.decim_rate = 2
    self.freq_band = [1., 40.]

