""" Selection of template events
  Input
    fpha: located phase file (with evid)
    fevt: use evid to get the event name
  Output
    template phase file
"""
from obspy import UTCDateTime

# i/o paths
fpha = 'input/eg_pal_hyp_full.pha'
fevt = 'input/eg_pal.evt'
ftemp = open('input/eg_pal.temp','w')
# selection criteria
ot_range = '20190704-20190707'
ot_range = [UTCDateTime(code) for code in ot_range.split('-')]
lat_range = [35.5,36.]
lon_range = [-117.8,-117.3]

print('get event_dict')
event_dict = {}
f=open(fevt); lines=f.readlines(); f.close()
for line in lines:
    evid, event_name = line.split(',')
    event_dict[evid] = event_name[:-1]

print('selecting template')
f=open(fpha); lines=f.readlines(); f.close()
for line in lines:
    codes = line.split(',')
    # if event line
    if len(codes[0])>=14:
        is_temp = True
        ot = UTCDateTime(codes[0])
        if not ot_range[0]<ot<ot_range[1]: is_temp = False
        lat, lon = [float(code) for code in codes[1:3]]
        if not lat_range[0]<lat<lat_range[1]: is_temp = False
        if not lon_range[0]<lon<lon_range[1]: is_temp = False
        evid = codes[-1][:-1]
        event_name = event_dict[evid]
        ot, lat, lon, dep, mag = codes[0:5]
        if is_temp: ftemp.write('{}_{},{},{},{},{},{}\n'.format(evid, event_name, ot, lat, lon, dep, mag))
    else:
        if is_temp: ftemp.write(line)

ftemp.close()
