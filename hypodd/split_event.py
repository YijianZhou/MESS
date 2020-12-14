""" Split event.dat if the number exceed the limit of hypoDD
"""
from obspy import UTCDateTime
import config

# i/o paths
fevent = 'input/event.dat'
# split params
cfg = config.Config()
time_ranges = cfg.time_ranges
start_id = cfg.evid_stride # starting evid for new detections

# read event.dat
event_list = []
f=open(fevent); lines=f.readlines(); f.close()
for line in lines:
    codes = line.split()
    datetime = UTCDateTime(codes[0]) + 1
    evid = int(codes[-1])
    event_list.append([evid, datetime, line])


for time_range in time_ranges:
    fout = open('input/event_%s.dat'%time_range, 'w')
    start_time, end_time = [UTCDateTime(code) for code in time_range.split('-')]
    for [evid, datetime, line] in event_list:
        if start_time<datetime<end_time or evid<start_id: 
            fout.write(line)
    fout.close()

