import unittest
import os
import uptide
import uptide.tidal_netcdf
import numpy
import datetime
try:
  import uptide.fes2012
except ImportError:
  pass

try:
  fes_data = os.environ['FES_DATA']
except:
  fes_data = None

test_ini_file = "tests/temp.ini"
trange = numpy.arange(0, 100)*100.
lat = 58.
lon = -6.

cnes_epoch = datetime.datetime(1950,1,1,0,0)

def create_ini_file(file_name, ini, constituents):
  """Create an ini file that is a subset of `ini` indicated by `constituents`"""
  f = file(file_name, 'w')
  for c in constituents:
    ci = ini['TIDE'][c]
    for key, value in ci.iteritems():
      f.write("TIDE_{0}_{1} = {2}\n".format(c,key,value))

  f.close()

if fes_data is not None:
  class TestFesUptideComparison(unittest.TestCase):
    def setUp(self):
      # read the ini file that contains all constituents available in FES2012
      # setting fes_data_path to "${FES_DATA}" means ${FES_DATA} is not expanded
      # but replaced with itself (i.e. unchanged) in the ini dictionary
      self.ini = uptide.tidal_netcdf.read_fes_ini_file(os.path.join(fes_data, "fes_2012.ini"), '${FES_DATA}')
      self.constituents = self.ini['TIDE'].keys()

    def test_individual_constituents(self):
      # fes adds additional long term corrections that we haven't implemented, so we only test short periods
      fes_long_period = ["MM", "MF", "MTM", "MSQM", "SSA"]
      # these constituents cause others to be infered (which we haven't implemented yet)
      fes_admittance_sources = ['K1', 'Q1', 'O1', 'M2', 'N2', 'K2', '2N2']
      fes_skip = fes_long_period+fes_admittance_sources
      for c in self.constituents:
        if c in fes_skip:
          continue

        create_ini_file(test_ini_file, self.ini, [c])
        # Compute tides using libfes:
        fes = uptide.fes2012.Fes(test_ini_file)
        h_fes = []; h_fes_short=[]
        for t in trange:
          h_fes.append(fes.get_tide(lat, lon, t/86400.))
          h_fes_short.append(fes.get_short_and_long_period_tide(lat, lon, t/86400.)[0])
        # Compute ties using the FES2012TidalInterpolator
        tide = uptide.Tides([c])
        tide.set_initial_time(cnes_epoch+datetime.timedelta(trange[0]))
        tnci = uptide.tidal_netcdf.FES2012TidalInterpolator(tide,
            test_ini_file, fes_data, ranges=[[360.+lon-1.0,360.+lon+1.0],[lat-1.0, lat+1.0]])
        h_uptide = []
        for t in trange:
          tnci.set_time(t)
          h_uptide.append(tnci.get_val((360.+lon, lat)))

        h_fes = numpy.array(h_fes)
        h_fes_short = numpy.array(h_fes_short)
        h_uptide = numpy.array(h_uptide)
        error = abs(h_fes_short-h_uptide).max()/abs(h_uptide+h_fes_short).max()
        if c == "L2":
          # not entirely sure what's wrong with L2
          self.assertLess(error, 0.3)
        else:
          self.assertLess(error, 0.007)


if __name__ == '__main__':
      unittest.main()
  
