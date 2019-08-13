import os
from obspy import UTCDateTime

# parallel params
start_date = UTCDateTime('20190704')
end_date   = UTCDateTime('20190713')
num_proc = 9
dt = (end_date - start_date) / num_proc

for i in range(num_proc):
    t0 = str((start_date + i*dt).date)
    t1 = str((start_date + (i+1)*dt).date)
    write_temp = 1 if i==0 else 0
    name_tag = '{}_{}'.format(t0,t1)
    start_id = 100000 * (i+1)
    os.system("python mk_pha.py \
        --time_range={},{} \
        --out_pha=./input2/pha_{} \
        --write_temp={} \
        --start_id={} &"\
        .format(t0, t1, name_tag, write_temp, start_id))

