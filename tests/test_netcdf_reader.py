import unittest
from uptide.netcdf_reader import NetCDFGrid, CoordinateError, NetCDFFile
import os
from numpy import arange, array, ones


# function used to fill the netcdf field, has to be linear
def f(lat, lon):
    return lat*10 + lon


test_file_name1 = 'tests/test_netcdf_reader1.nc'
test_file_name2 = 'tests/test_netcdf_reader2.nc'


class TestNetCDFGrid(unittest.TestCase):
    """Tests the uptide.netcdf_reader.NetCDFGrid class"""
    def setUp(self):
        # it seems that many scipy installations are broken for
        # netcdf writing - therefore simply committing the
        # test files instead of writing them out on the fly here
        return
        zval = array([[f(lat, lon) for lon in arange(10.0)]
                      for lat in arange(10.0)])
        nc = NetCDFFile(test_file_name1, 'w')
        nc.createDimension('lat', 10)
        nc.createDimension('lon', 10)
        nc.createVariable('latitude', 'float64', ('lat', ))
        nc.createVariable('longitude', 'float64', ('lon', ))
        nc.variables['latitude'][:] = arange(10.0)
        nc.variables['longitude'][:] = arange(10.0)
        nc.createVariable('z', 'float64', ('lat', 'lon'))
        nc.variables['z'][:, :] = zval
        nc.createVariable('mask', 'float64', ('lat', 'lon'))
        mask = ones((10, 10), dtype='float64')
        mask[0:2, :] = 0.0
        nc.variables['mask'][:, :] = mask
        nc.createVariable('transposed_mask', 'float64', ('lon', 'lat'))
        nc.variables['transposed_mask'][:, :] = mask.T
        nc.close()
        # same thing but without the coordinate fields and mask
        nc = NetCDFFile(test_file_name2, 'w')
        nc.createDimension('lat', 10)
        nc.createDimension('lon', 10)
        nc.createVariable('z', 'float64', ('lat', 'lon'))
        nc.variables['z'][:, :] = zval
        nc.close()

    def tearDown(self):
        # don't remove them either (see above)
        return
        os.remove(test_file_name1)
        os.remove(test_file_name2)
        pass

    def _test_interpolator(self, interpolator, mask, ranges, coordinate_perm):
        # first the tests common to all combinations
        # point that is always inside:
        xy = [[4.33, 5.2][i] for i in coordinate_perm]
        self.assertEqual(interpolator.get_val(xy), f(4.33, 5.2))
        # point outside the domain, should raise exception:
        xy = [[-4.95, 8.3][i] for i in coordinate_perm]
        self.assertRaises(CoordinateError, interpolator.get_val, xy)

        if mask is not None:
            # point within sea, should work as before:
            xy = [[4.33, 5.2][i] for i in coordinate_perm]
            self.assertAlmostEqual(interpolator.get_val(xy), f(4.33, 5.2))
            # point between row of land and of sea points, should interpolate from nearest sea row:
            xy = [[1.2, 8.3][i] for i in coordinate_perm]
            self.assertAlmostEqual(interpolator.get_val(xy), f(2.0, 8.3))
            # point inside the first two land rows, should raise exception
            xy = [[0.95, 8.3][i] for i in coordinate_perm]
            self.assertRaises(CoordinateError, interpolator.get_val, xy)
        if ranges:
            # test within the range
            xy = [[2.9, 7.0][i] for i in coordinate_perm]
            self.assertAlmostEqual(interpolator.get_val(xy), f(2.9, 7.))
            # tests outside the range, should raise exception
            xy = [[3.2, 0.9][i] for i in coordinate_perm]
            self.assertRaises(CoordinateError, interpolator.get_val, xy)
            xy = [[5.9, 9.0][i] for i in coordinate_perm]
            self.assertRaises(CoordinateError, interpolator.get_val, xy)

    # test a specific combination of mask (None, 'mask', 'transposed_mask' or 'mask_from_fill_value'),
    # ranges (True or False) and coordinate_perm (lat, lon) or (lon, lat)
    def _test_combination(self, mask, ranges, coordinate_perm):
        kwargs = {}
        if mask == 'mask_from_fill_value':
            kwargs['mask_from_fill_value'] = ('mask', 0.0)
        elif mask is not None:
            kwargs['mask'] = mask
        if ranges:
            kwargs['ranges'] = [((0., 4.), (2., 8.))[i] for i in coordinate_perm]
        # load the netcdf created in setup()
        if coordinate_perm == (0, 1):
            grid = NetCDFGrid(test_file_name1, ('lat', 'lon'), ('latitude', 'longitude'), **kwargs)
        else:
            grid = NetCDFGrid(test_file_name1, ('lon', 'lat'), ('longitude', 'latitude'), **kwargs)
        self._test_interpolator(grid.get_interpolator('z'), mask, ranges, coordinate_perm)

        # now try the same for the case where the field values are stored in a separate file
        grid2 = NetCDFGrid(test_file_name2, grid)
        self._test_interpolator(grid2.get_interpolator('z'), mask, ranges, coordinate_perm)

    # test all combinations of mask and ranges options
    # also try out coordinate permutations lat, lon and lon, lat (the read nc file is lat, lon in both cases)
    def test_all_combinations(self):
        for mask in (None, 'mask', 'transposed_mask', 'mask_from_fill_value'):
            for ranges in (False, True):
                for coordinate_perm in ((0, 1), (1, 0)):
                    self._test_combination(mask, ranges, coordinate_perm)

    def test_get_values(self):
        grid = NetCDFGrid(test_file_name1, ('lon', 'lat'), ('longitude', 'latitude'), ranges=((2., 8.), (0., 4.)))
        val = grid.get_values('z')
        self.assertEqual(val.shape, tuple(grid.shape))
        lon0, lat0 = grid.origin
        self.assertAlmostEqual(val[1, 2], f(lat0 + 2., lon0 + 1.))


if __name__ == '__main__':
    unittest.main()
