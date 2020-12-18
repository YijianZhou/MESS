import os, shutil

# i/o paths
gpu_idx = '1'
mess_dir = '/home/zhouyj/software/MESS'
data_dir = '/data2/Example_data'
time_range = '20190705-20190710'
sta_file = 'input/example_pad.sta'
temp_root = '/data3/bigdata/zhouyj/Example_templates'
temp_pha = 'input/example_pad.temp'
out_root = 'output/example'
out_ctlg = os.path.join(out_root, 'catalog_{}.dat'.format(time_range))
out_pha = os.path.join(out_root, 'phase_{}.dat'.format(time_range))

# run MESS (gpu ver)
shutil.copyfile('config_example.py', os.path.join(mess_dir, 'config.py'))
os.system("python {}/run_mess_gpu.py --gpu_idx={} \
    --data_dir={}  --time_range={} --sta_file={} \
    --temp_root={} --temp_pha={} \
    --out_ctlg={} --out_pha={}"\
    .format(mess_dir, gpu_idx, data_dir, time_range, sta_file, temp_root, temp_pha, out_ctlg, out_pha))

