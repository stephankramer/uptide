import unittest
import uptide
from uptide.netcdf_reader import NetCDFInterpolator, CoordinateError
import itertools
import os
from scipy.io.netcdf import NetCDFFile
from numpy import arange, array, ones

# function used to fill the netcdf field, has to be linear
def f(lat, lon):
  return lat*10 + lon

class TestNetCDFInterpolator(unittest.TestCase):
  """Tests the uptide.netcdf.NetCDFInterpolator class"""
  def setUp(self):
    zval = array(
        [[f(lat,lon) for lon in arange(10.0)]
                     for lat in arange(10.0)])
    nc = NetCDFFile('test.nc', 'w')
    nc.createDimension('lat', 10)
    nc.createDimension('lon', 10)
    nc.createVariable('latitude', 'float64', ('lat',))
    nc.createVariable('longitude', 'float64', ('lon',))
    nc.variables['latitude'][:] = arange(10.0)
    nc.variables['longitude'][:] = arange(10.0)
    nc.createVariable('z', 'float64', ('lat','lon'))
    nc.variables['z'][:,:] = zval
    nc.createVariable('mask', 'float64', ('lat','lon'))
    mask = ones((10,10),dtype='float64')
    mask[0:2,:] = 0.0
    nc.variables['mask'][:,:] = mask
    nc.close()
    # same thing but without the coordinate fields and mask
    nc = NetCDFFile('test2.nc', 'w')
    nc.createDimension('lat', 10)
    nc.createDimension('lon', 10)
    nc.createVariable('z', 'float64', ('lat','lon'))
    nc.variables['z'][:,:] = zval
    nc.close()

  def tearDown(self):
    os.remove('test.nc')
    os.remove('test2.nc')

  def _test_prepared_nci(self, nci, perm):
    # first the tests common to all permutations
    # point that is always inside:
    self.assertEqual(nci.get_val((4.33,5.2)), f(4.33,5.2))
    # point outside the domain, should raise exception:
    self.assertRaises(CoordinateError, nci.get_val, (-4.95,8.3))

    if 'mask' in perm or 'mask_from_fill_value' in perm:
      # point within sea, should work as before:
      self.assertAlmostEqual(nci.get_val((4.33,5.2)), f(4.33,5.2))
      # point between row of land and of sea points, should interpolate from nearest sea row:
      self.assertAlmostEqual(nci.get_val((1.2,8.3)), f(2.0,8.3))
      # point inside the first two land rows, should raise exception
      self.assertRaises(CoordinateError, nci.get_val, (0.95,8.3))
    if 'ranges' in perm:
      # test within the range
      self.assertAlmostEqual(nci.get_val((2.9,7.0)), f(2.9,7.))
      # tests outside the range, should raise exception
      self.assertRaises(CoordinateError, nci.get_val, (3.2,0.9))
      self.assertRaises(CoordinateError, nci.get_val, (5.9,9.))

  # test a specific permutation of the calling sequence set_field, set_mask, set_ranges
  def _test_permutation(self, perm):
    # load the netcdf created in setup()
    nci = NetCDFInterpolator('test.nc', ('lat', 'lon'), ('latitude', 'longitude'))
    # call the methods in the order given by perm
    for x in perm:
      if x=='field':
        nci.set_field('z')
      elif x=='mask':
        nci.set_mask('mask')
      elif x=='mask_from_fill_value':
        nci.set_mask_from_fill_value('mask', 0.0)
      elif x=='ranges':
        nci.set_ranges(((0.,4.),(2.,8.)))
      else:
        raise Exception("Unknown method")

    # now perform all tests
    self._test_prepared_nci(nci, perm)

    # now try the same for the case where the field values are stored in a separate file
    nci2 = NetCDFInterpolator('test2.nc', nci)
    nci2.set_field('z')
    self._test_prepared_nci(nci2, perm)


  # test all permutations of the calling sequence set_field, set_mask, set_ranges
  # including all permutations that only call 1 or 2 of these methods
  # set_field should always be called
  def test_all_permutations(self):
    for n in range(1,4):
      for perm in itertools.permutations(['field','mask','ranges'], n):
        if not 'field' in perm:
          continue
        self._test_permutation(perm)

  def test_all_permutations_with_fill_value(self):
    for n in range(1,4):
      for perm in itertools.permutations(['field','mask_from_fill_value','ranges'], n):
        if not 'field' in perm:
          continue
        self._test_permutation(perm)
