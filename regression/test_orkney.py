from __future__ import print_function
import uptide
import uptide.tidal_netcdf
import os.path
import datetime
import numpy
import utm
import sys

utm_zone = 30
utm_band = 'V'
ranges = ((-4.0,0.0), (58.0,61.0))
otps_data_dir = '/home/skramer/git/OpenTidalFarm/examples/orkney_resource_assessment'
grid_file_name = os.path.join(otps_data_dir, 'gridES2008.nc')
data_file_name = os.path.join(otps_data_dir, 'hf.ES2008.nc')
initial_time = datetime.datetime(2001, 9, 18, 0)
constituents = ['Q1', 'O1', 'P1', 'K1', 'N2', 'M2', 'S2', 'K2']

reference_file_name = 'orkney.dat'

times = (0., 100.)
boundary_xy = numpy.loadtxt('orkney.xy')

tide = uptide.Tides(constituents)
tide.set_initial_time(initial_time)
tnci = uptide.tidal_netcdf.OTPSncTidalInterpolator(tide,
            grid_file_name, data_file_name, ranges)

vals_t = []
for t in times:
    tnci.set_time(t)
    vals = []
    for xy in boundary_xy:
        latlon = utm.to_latlon(xy[0], xy[1], utm_zone, utm_band)
        vals.append(tnci.get_val((latlon[1], latlon[0]), allow_extrapolation=True))

    vals_t.append(vals)

vals_t = numpy.array(vals_t).T
# uncomment to update reference:
#numpy.savetxt(reference_file_name, vals_t)

reference = numpy.loadtxt(reference_file_name)
err =  numpy.abs(reference-vals_t).max()
print("err = ", err)
if err>1e-3:
    print("Error too large")
    sys.exit(1)
else:
    print("Success!")
