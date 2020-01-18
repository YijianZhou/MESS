import os, shutil

# i/o paths
msms_dir = '/home/zhouyj/software/MSMS'
data_dir = '/data2/ZSY_SAC/*/*'
date_rng = '20160901-20190201'
out_root = './output/zsy'
temp_root = os.path.join(out_root, 'Templates')
temp_pha = os.path.join(out_root, 'zsy_201609-201901_hyp.pha')
out_ctlg = os.path.join(out_root, 'phase_{}.dat'.format(date_rng))
out_pha = os.path.join(out_root, 'catalog_{}.dat'.format(date_rng))

# run MSMS cut temp
shutil.copyfile('config.py', os.path.join(msms_dir, 'config.py'))
os.system("python {}/run_msms_gpu.py \
    --data_dir={} \
    --date_range={} \
    --temp_root={} \
    --temp_pha={} \
    --out_ctlg={} \
    --out_pha={}"\
    .format(msms_dir, data_dir, date_rng, temp_root, temp_pha, out_ctlg, out_pha))

