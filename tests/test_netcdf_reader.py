import unittest
import pytest
import netCDF4
import numpy as np
from uptide.netcdf_reader import NetCDFGrid, CoordinateError, NetCDFFile, NetCDFInterpolatorError
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


# tests for grids with periodic longitude: field values are g(lon, lat) in the grid points
def g(lon, lat):
    return np.cos(np.radians(lon)) + 0.1*lat


def write_lonlat_file(path, lons, lats, transposed=False, land_lon=None):
    """Write field z = g(lon, lat) on a lon, lat grid. If transposed, store as (lat, lon).
    If land_lon is given, the column with that longitude is masked out in the 'mask' field."""
    ds = netCDF4.Dataset(path, 'w')
    ds.createDimension('nx', len(lons))
    ds.createDimension('ny', len(lats))
    ds.createVariable('lon', 'float64', ('nx',))[:] = lons
    ds.createVariable('lat', 'float64', ('ny',))[:] = lats
    lon2d, lat2d = np.meshgrid(lons, lats, indexing='ij')
    mask = np.ones_like(lon2d)
    if land_lon is not None:
        mask[np.isclose(lon2d, land_lon)] = 0.
    dims = ('ny', 'nx') if transposed else ('nx', 'ny')
    for name, val in [('z', g(lon2d, lat2d)), ('mask', mask)]:
        ds.createVariable(name, 'float64', dims)[:] = val.T if transposed else val
    ds.close()
    return str(path)


def expected_val(lon, lat, lons):
    """Expected interpolated value at lon, lat, interpolating linearly between the neighbouring
    grid longitudes lons (taken modulo 360)."""
    dlon = lons[1] - lons[0]
    lon = lons[0] + (lon - lons[0]) % 360.
    l0 = lons[0] + dlon*np.floor((lon - lons[0])/dlon)
    alpha = (lon - l0)/dlon
    return (1-alpha)*g(l0, lat) + alpha*g(l0+dlon, lat)


def open_lonlat_grid(path, transposed_dims, **kwargs):
    if transposed_dims:
        if 'ranges' in kwargs:
            kwargs['ranges'] = kwargs['ranges'][::-1]
        grid = NetCDFGrid(path, ('ny', 'nx'), ('lat', 'lon'), coordinate_periods=(None, 360.), **kwargs)
        interpolator = grid.get_interpolator('z')
        return grid, lambda lon, lat, *args: interpolator.get_val((lat, lon), *args)
    else:
        grid = NetCDFGrid(path, ('nx', 'ny'), ('lon', 'lat'), coordinate_periods=(360., None), **kwargs)
        interpolator = grid.get_interpolator('z')
        return grid, lambda lon, lat, *args: interpolator.get_val((lon, lat), *args)


lats = np.array([-10., 0., 10.])
# grid covering exactly one period (like TPXO9)
periodic_lons = np.arange(1., 13.)*30.
# grid covering one period with first and last longitude repeated
repeated_lons = np.arange(0., 13.)*30.


@pytest.mark.parametrize('lons', [periodic_lons, repeated_lons], ids=['periodic', 'repeated'])
@pytest.mark.parametrize('transposed', [False, True])
@pytest.mark.parametrize('transposed_dims', [False, True])
@pytest.mark.parametrize('ranges', [None, ((-40., 40.), (-5., 5.)), ((350., 370.), (-5., 5.)),
                                    ((-180., 180.), (-5., 5.)), ((100., 200.), (-5., 5.))])
def test_periodic_grid(tmp_path, lons, transposed, transposed_dims, ranges):
    path = write_lonlat_file(tmp_path / 'grid.nc', lons, lats, transposed)
    kwargs = {} if ranges is None else {'ranges': ranges}
    grid, get_val = open_lonlat_grid(path, transposed_dims, **kwargs)
    if ranges is None or ranges[0][1] - ranges[0][0] >= 360.:
        test_lons = [0., 15., -15., 345., 375., 200., -705., 360.]
    else:
        test_lons = [x + offset for x in np.linspace(*ranges[0], 7) for offset in (-360., 0., 360.)]
    for lon in test_lons:
        assert np.isclose(get_val(lon, 2.5), expected_val(lon, 2.5, lons)), lon
    # same values when reading the values from a separate file
    grid2 = NetCDFGrid(path, grid)
    get_val2 = grid2.get_interpolator('z').get_val
    lon = test_lons[1]
    x = (2.5, lon) if transposed_dims else (lon, 2.5)
    assert np.isclose(get_val2(x), expected_val(lon, 2.5, lons))


@pytest.mark.parametrize('transposed', [False, True])
def test_periodic_grid_mask(tmp_path, transposed):
    # mask out the last column (lon=360), so that a point between 360 and 30 only interpolates from 30
    path = write_lonlat_file(tmp_path / 'grid.nc', periodic_lons, lats, transposed, land_lon=360.)
    extrapolated = []
    for kwargs in [{'mask': 'mask'}, {'mask': 'mask', 'ranges': ((-40., 40.), (-5., 5.))}]:
        grid, get_val = open_lonlat_grid(path, False, **kwargs)
        assert np.isclose(get_val(15., 2.5), g(30., 2.5))
        assert np.isclose(get_val(-15., 2.5), g(330., 2.5))
        with pytest.raises(CoordinateError):
            get_val(0., 2.5)
        # extrapolation uses sea points on both sides of the seam, so should not depend on ranges
        extrapolated.append(get_val(0., 2.5, True))
    assert np.isclose(*extrapolated)


def test_periodic_grid_values(tmp_path):
    path = write_lonlat_file(tmp_path / 'grid.nc', periodic_lons, lats, transposed=True)
    grid, _ = open_lonlat_grid(path, False, ranges=((-40., 40.), (-5., 5.)))
    val = grid.get_values('z')
    assert val.shape == tuple(grid.shape)
    lon = grid.origin[0] + grid.delta[0]*np.arange(grid.shape[0])
    lat = grid.origin[1] + grid.delta[1]*np.arange(grid.shape[1])
    np.testing.assert_allclose(val, g(*np.meshgrid(lon, lat, indexing='ij')), atol=1e-12)


regional_lons = np.arange(-60., 61., 10.)


@pytest.mark.parametrize('transposed_dims', [False, True])
def test_regional_grid(tmp_path, transposed_dims):
    path = write_lonlat_file(tmp_path / 'grid.nc', regional_lons, lats)
    for ranges in [None, ((330., 350.), (-5., 5.)), ((-30., -10.), (-5., 5.))]:
        kwargs = {} if ranges is None else {'ranges': ranges}
        _, get_val = open_lonlat_grid(path, transposed_dims, **kwargs)
        for lon in [-25., 335., 695., -385.]:
            assert np.isclose(get_val(lon, 2.5), expected_val(lon, 2.5, regional_lons)), lon
        for lon in [180., -65., 65., 295.]:
            with pytest.raises(CoordinateError):
                get_val(lon, 2.5)
    # range only partly overlapping the west edge of the grid
    _, get_val = open_lonlat_grid(path, transposed_dims, ranges=((-70., -50.), (-5., 5.)))
    assert np.isclose(get_val(-55., 2.5), expected_val(-55., 2.5, regional_lons))
    with pytest.raises(NetCDFInterpolatorError):
        open_lonlat_grid(path, transposed_dims, ranges=((100., 200.), (-5., 5.)))


def test_regional_grid_across_dateline(tmp_path):
    lons = np.arange(150., 211., 10.)
    path = write_lonlat_file(tmp_path / 'grid.nc', lons, lats)
    # ranges as computed from min and max of -180<lon<180 for a domain crossing the dateline
    for ranges in [((-180., 180.), (-5., 5.)), ((-179., 179.), (-5., 5.))]:
        grid, get_val = open_lonlat_grid(path, False, ranges=ranges)
        assert grid.shape[0] == len(lons)
        for lon in [155., 175., -175., -155.]:
            assert np.isclose(get_val(lon, 2.5), expected_val(lon, 2.5, lons)), lon


def test_nearly_periodic_grid_warns(tmp_path):
    lons = np.arange(1., 13.)*360./11.5
    path = write_lonlat_file(tmp_path / 'grid.nc', lons, lats)
    with pytest.warns(UserWarning, match='nearly'):
        _, get_val = open_lonlat_grid(path, False)
    with pytest.raises(CoordinateError):
        get_val(lons[0] - 1., 2.5)


if __name__ == '__main__':
    unittest.main()
