import os
from obspy import UTCDateTime

# parallel params
date_rng = '20180325-20190201'
num_proc = 20
start_date, end_date = [UTCDateTime(date) for date in date_rng.split('-')]
dt = (end_date - start_date) / num_proc

for i in range(num_proc):
    t0 = ''.join(str((start_date + i*dt).date).split('-'))
    t1 = ''.join(str((start_date + (i+1)*dt).date).split('-'))
    date_rng = '{}-{}'.format(t0,t1)
    out_dt = 'input/dt_{}'.format(date_rng)
    out_event = 'input/event_{}'.format(date_rng)
    start_id = 100000 * (i+1)
    os.system("python mk_dt.py \
        --date_range={} --start_id={} \
        --out_dt={} --out_event={} &" \
        .format(date_rng, start_id, out_dt, out_event))

