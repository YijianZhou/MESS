t0='20190707'
t1='20190710'
net='rc'
python run_mft_gpu.py --time_range=$t0,$t1 \
--data_dir=/data2/Ridgecrest/*/* \
--out_pha=./output/$net/phase_${t0}_${t1}_aug.dat \
--out_ctlg=./output/$net/catalog_${t0}_${t1}_aug.dat \
--temp_root=./output/rc/Templates \
--temp_pha=./output/rc/phase_rc1.dat

