import os
from obspy import UTCDateTime

# parallel params
start_date = UTCDateTime('20170927')
end_date   = UTCDateTime('20171219')
num_proc = 20
dt = (end_date - start_date) / num_proc

for i in range(num_proc):
    t0 = str((start_date + i*dt).date)
    t1 = str((start_date + (i+1)*dt).date)
    name_tag = '{}_{}'.format(t0,t1)
    start_id = 100000 * (i+1)
    os.system("python mk_dt.py \
        --time_range={},{} \
        --out_dt=./input/dt_{} \
        --out_event=./input/event_{} \
        --start_id={} &"\
        .format(t0, t1, name_tag, name_tag, start_id))

