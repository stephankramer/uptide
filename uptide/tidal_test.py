import tidal_netcdf
import tidal_composer
from datetime import datetime, timedelta
from numpy import arange, loadtxt
from math import cos, pi
from pylab import plot, show, figure
import os
import time

constituents= ('M2','S2','N2','K2','K1','O1','P1','Q1')
lat = 58.708
lon = -3.287
dt0 = datetime(2003,3,28,0,0,0)
trange = arange(0.,24.*32., 0.25)*3600.

tide = tidal_composer.Tides(constituents)
tide.set_initial_time(dt0)

tn = tidal_netcdf.AMCGTidalInterpolator(tide,
        '/home/skramer/bld/OTPSnc/ES2008/DATA/ap.ES2008.nc',
        ranges=((58.0,61.0),(-4.0,0.0)))

etas = []
for t in trange:
  tn.set_time(t)
  eta = tn.get_val((lat,lon))
  etas.append(eta)

figure()
plot(trange, etas)

tnci=tidal_netcdf.OTPSncTidalInterpolator(tide,
         '/home/skramer/bld/OTPSnc/ES2008/DATA/gridES2008.nc',
         '/home/skramer/bld/OTPSnc/ES2008/DATA/hf.ES2008.nc',
        ranges=((-4.0,0.0),(58.0,61.0)))
etas2 = []
for t in trange:
  tnci.set_time(t)
  eta = tnci.get_val((lon,lat))
  etas2.append(eta)

plot(trange, etas2)

tnci=tidal_netcdf.FESTidalInterpolator(tide,
         '/home/skramer/data/tide.fes2004.nc',
        ranges=((58.0,61.0),(356.,360.0)))
etas3 = []
for t in trange:
  tnci.set_time(t)
  eta = tnci.get_val((lat,360.+lon))
  etas3.append(eta)

plot(trange, etas3)

f = file('latlontime.dat', 'w')
for t in trange:
  dt = dt0 + timedelta(seconds=t)
  f.write('%f %f '%(lat,lon) + dt.strftime('%Y %m %d %H %M %S')+'\n')
f.close()

os.system('/home/skramer/bld/OTPSnc/predict_tide < setup.inp')
otps_eta=loadtxt('sample.out',skiprows=6,usecols=(4,))

plot(trange, otps_eta)
show()
