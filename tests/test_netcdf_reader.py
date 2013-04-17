import unittest
import uptide
import uptide.netcdf_reader as un
from scipy.io.netcdf import NetCDFFile

class TestNetCDFInterpolator(unittest.TestCase):
  """Tests the uptide.netcdf.NetCDFInterpolator class"""
  def setUp(self):
    nc = NetCDFFile('test.nc', 'w')
