import uptide
import uptide.tidal_netcdf
import datetime
from numpy import *

constituents = ['M2', 'S2']
tide = uptide.Tides(constituents)
tide.set_initial_time(datetime.datetime(1950,1,1,0,0))
tnci = uptide.tidal_netcdf.FES2012TidalInterpolator(tide, '/data/packages/fes2012/data/fes_2012.ini')

data = loadtxt('/data/packages/fes2012/examples/fes_slev.txt', skiprows=1)
for t in data[:,0]:
  tnci.set_time(datetime.timedelta(t).total_seconds())
  print tnci.get_val((360.+-7.688,59.195))
