import os, shutil

# i/o paths
msms_dir = '/home/zhouyj/software/MSMS'
data_dir = '/data2/ZSY_SAC/*/*'
out_root = 'input/zsy'
temp_pha = os.path.join(out_root, 'zsy_201609-201901_hyp.pha')

# run MSMS cut temp
shutil.copyfile('config.py', os.path.join(msms_dir, 'config.py'))
os.system("python {}/cut_temp.py \
    --data_dir={} --temp_pha={} --out_root={}" \
    .format(msms_dir, data_dir, temp_pha, out_root))

