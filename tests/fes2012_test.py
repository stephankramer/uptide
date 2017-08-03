from __future__ import print_function
import uptide
import uptide.tidal_netcdf
import datetime
from numpy import *
import os
import pytest

pytest.importorskip('uptide.fes2012')

fes_dir='/data/packages/fes2012'
fes_dir='/home/skramer/data/fes2012'
data = loadtxt(fes_dir+'/examples/fes_slev.txt', skiprows=1)
cnes_epoch = datetime.datetime(1950,1,1,0,0)
lon=-7.688; lat=59.195;
constituents = ['M2', 'S2', 'N2', 'K2',
  'Q1', 'O1', 'K1', 'P1',
  'M4', 'S1', '2N2',
  'MTM', 'MM', 'MF', 'MSQM']
#constituents = ['MSQM',]
tide = uptide.Tides(constituents)
tide.set_initial_time(cnes_epoch+datetime.timedelta(data[0,0]))
#tide.f[:]=1.0
#tide.u[:]=0.0
tnci = uptide.tidal_netcdf.FES2012TidalInterpolator(tide, fes_dir+'/data/fes.ini')

ut = []
for t in data[:,0]:
  tnci.set_time((t-data[0,0])*86400.)
  ut.append(tnci.get_val((360.+lon,lat)))

f = file(fes_dir+'/data/test.ini', 'w')
for c in constituents:
  f.write("""TIDE_{0}_AMPLITUDE         = amplitude
TIDE_{0}_FILE              = ${{FES_DATA}}/{0}_tide.nc
TIDE_{0}_LATITUDE          = latitude
TIDE_{0}_LONGITUDE         = longitude
TIDE_{0}_PHASE             = phase
TIDE_{0}_PHASE_LAG         = 1
""".format(c))

f.close()


os.system('cd {}/examples_hacked/; ./fes_slev > test.txt'.format(fes_dir))
fes = loadtxt(fes_dir+'/examples_hacked/test.txt', skiprows=1, usecols=(7,))
err = abs(fes-ut).max()
print(constituents, ':',  err, err/max(abs(fes)))
