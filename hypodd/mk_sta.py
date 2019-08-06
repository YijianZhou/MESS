sta_file = '/home/zhouyj/Desktop/California/preprocess/station.dat'
f=open(sta_file); lines=f.readlines(); f.close()
out = open('station.dat','w')

for line in lines:
    net, sta, chn, lon, lat, ele = line.split(',')
    lon = float(lon)
    lat = float(lat)
    out.write('{} {} {}\n'.format(sta, lat, lon))
out.close()
