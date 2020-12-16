""" Selection of template events
  Input
    hypoInverse catalog
    fevt: use evid to get the event name
  Output
    template phase file
"""
# i/o paths
fpha = 'input/example_pad_hyp_all.pha' # --> event loc
fevt = 'input/example_pad.evt' # --> event name
ftemp = open('input/example_pad.temp','w')
# selection criteria
lat_rng = [35.5,36.]
lon_rng = [-117.8,-117.3]

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
        lat, lon = [float(code) for code in codes[1:3]]
        if not lat_rng[0]<lat<lat_rng[1]: is_temp = False
        if not lon_rng[0]<lon<lon_rng[1]: is_temp = False
        evid = codes[-1][:-1]
        event_name = event_dict[evid]
        ot, lat, lon, dep, mag = codes[0:5]
        if is_temp: ftemp.write('{}_{},{},{},{},{},{}\n'.format(evid, event_name, ot, lat, lon, dep, mag))
    else:
        if is_temp: ftemp.write(line)

ftemp.close()
