import config

# i/o paths
cfg = config.Config()
fsta = cfg.fsta
f=open(fsta); lines=f.readlines(); f.close()
out = open('input/station.dat','w')

for line in lines:
    net, sta, chn, lon, lat, ele = line.split(',')
    lon = float(lon)
    lat = float(lat)
    out.write('{} {} {}\n'.format(sta, lat, lon))
out.close()
