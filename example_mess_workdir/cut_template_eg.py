""" Cut template data
"""
import os, shutil

# i/o paths
mess_dir = '/home/zhouyj/software/MESS'
data_dir = '/data1/Example_data'
out_root = 'output/Example_templates'
temp_pha = 'input/eg_pal.temp'

shutil.copyfile('config_eg.py', os.path.join(mess_dir, 'config.py'))
os.system("python {}/cut_template_torch.py \
    --data_dir={} --temp_pha={} --out_root={}"\
    .format(mess_dir, data_dir, temp_pha, out_root))

