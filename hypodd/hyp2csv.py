def get_ot(info):
    yr  = info[10]
    mon = info[11]
    day = info[12]
    hr  = info[13]
    mn  = info[14]
    sec = info[15]
    return '{}{:0>2}{:0>2}{:0>2}{:0>2}{:0>6}'.format(yr, mon, day, hr, mn, sec)

hypodd_out = 'output/hypoDD.reloc'
f=open(hypodd_out); lines=f.readlines(); f.close()
out_csv = open('output/mft_reloc.csv','w')
for line in lines:
    info = line.split()
    evid = float(info[0])
    if evid<10000: continue #TODO
    lat = info[1]
    lon = info[2]
    dep = info[3]
    mag = info[16]
    ot = get_ot(info)
    out_csv.write('{},{},{},{},{}\n'.format(ot, lat, lon, dep, mag))
out_csv.close()
