import os, shutil
from obspy import UTCDateTime

# i/o paths
gpu_idx = "0"
msms_dir = '/home/zhouyj/software/MSMS'
data_dir = '/data2/ZSY_SAC/*/*'
date_rng = '20160901-20190201'
num_days = 60 # per file
out_root = 'output/zsy'
temp_root = 'input/zsy/Templates'
temp_pha = 'input/zsy/zsy_201609-201901_hyp.pha'

# run MSMS cut temp
shutil.copyfile('config.py', os.path.join(msms_dir, 'config.py'))
start_date, end_date = [UTCDateTime(date) for date in date_rng.split('-')]
num_files = int((end_date - start_date) / (num_days * 86400)) + 1
for i in range(num_files):
    t0 = ''.join(str((start_date + i*num_days*86400).date).split('-'))
    t1 = ''.join(str((start_date + (i+1)*num_days*86400).date).split('-'))
    date_rng = '{}-{}'.format(t0, t1)
    out_pha = '{}/phase_{}.dat'.format(out_root, date_rng)
    out_ctlg = '{}/catalog_{}.dat'.format(out_root, date_rng)
    os.system("python {}/run_msms_gpu.py --gpu_idx={} \
        --data_dir={} --date_range={} \
        --temp_root={} --temp_pha={} \
        --out_ctlg={} --out_pha={}"\
        .format(msms_dir, gpu_idx, data_dir, date_rng, temp_root, temp_pha, out_ctlg, out_pha))

