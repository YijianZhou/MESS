""" Make station input file for HypoDD
"""
import os
import config

# i/o paths
cfg = config.Config()
fout = open('input/station.dat', 'w')

f=open(cfg.fsta); lines=f.readlines(); f.close()
for line in lines:
    net_sta, lat, lon, ele = line.split(',')[0:4]
    net, sta = net_sta.split('.')
    lon = float(lon)
    lat = float(lat)
    fout.write('{} {} {}\n'.format(sta, lat, lon))
fout.close()
