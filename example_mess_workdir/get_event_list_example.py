""" Get name list from phase file
  Input
    fpha: phase file of the original detection
  Output
    fevt: lines of 'evid,event_name'
"""
from obspy import UTCDateTime

# i/o paths
fpha = 'input/example_pad.pha'
fevt = open('input/example_pad.evt','w')

def dtime2str(dtime):
    date = ''.join(str(dtime).split('T')[0].split('-'))
    time = ''.join(str(dtime).split('T')[1].split(':'))[0:9]
    return date + time

# read phase file
evid = 0
f=open(fpha); lines=f.readlines(); f.close()
for line in lines:
    codes = line.split(',')
    if len(codes)!=5: continue
    event_name = dtime2str(UTCDateTime(codes[0]))
    fevt.write('{},{}\n'.format(evid, event_name))
    evid += 1

fevt.close()
