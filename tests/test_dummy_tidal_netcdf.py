import pytest
import uptide
import netCDF4
import os
import numpy as np
import datetime
from uptide.tidal_netcdf import TPXOncTidalComponentInterpolator

constituents = ('M2', 'S2', 'N2', 'K2', 'K1', 'O1', 'P1', 'Q1')


@pytest.fixture
def dummy_tpxo_grid_file(tmp_path):
    ds = netCDF4.Dataset(tmp_path / 'grid_tpxo9.nc', 'w', format='NETCDF3_CLASSIC')
    ds.createDimension('nx', 2)
    ds.createDimension('ny', 2)
    mz = ds.createVariable('mz', np.int32, ('nx', 'ny'))
    mz[:] = 1
    mu = ds.createVariable('mu', np.int32, ('nx', 'ny'))
    mu[:] = 1
    mv = ds.createVariable('mv', np.int32, ('nx', 'ny'))
    mv[:] = 1
    lon = [[0, 0], [360, 360]]
    lat = [[-90, 90], [-90, 90]]
    for c in 'zuv':
        lon_c = ds.createVariable(f'lon_{c}', np.int32, ('nx', 'ny'))
        lon_c[:] = lon
        lat_c = ds.createVariable(f'lat_{c}', np.int32, ('nx', 'ny'))
        lat_c[:] = lat
    filepath = ds.filepath()
    ds.close()
    yield filepath
    os.remove(filepath)


@pytest.fixture
def dummy_tpxo_elev_file(tmp_path):
    ds = netCDF4.Dataset(tmp_path / 'h_tpxo9.nc', 'w', format='NETCDF3_CLASSIC')
    ds.createDimension('nx', 2)
    ds.createDimension('ny', 2)
    ds.createDimension('nc', len(constituents))
    nct = ds.createDimension('nct', 4)
    con = ds.createVariable('con', 'c', ('nc', 'nct'))
    con[:] = [x.ljust(nct.size) for x in constituents]
    hRe = ds.createVariable('hRe', np.float64, ('nc', 'nx', 'ny'))
    # give M2 amplitude of 1, S2 amplitude of 2, etc.
    hRe[:] = (np.arange(len(constituents))[:, np.newaxis, np.newaxis] + 1.)
    hIm = ds.createVariable('hIm', np.float64, ('nc', 'nx', 'ny'))
    hIm[:] = np.zeros((len(constituents), 2, 2))
    filepath = ds.filepath()
    ds.close()
    yield filepath
    os.remove(filepath)


def test_dummy_tpxo(dummy_tpxo_grid_file, dummy_tpxo_elev_file):
    tide = uptide.Tides(['M2'])
    dt0 = datetime.datetime(2003, 3, 28, 0, 0, 0)
    tide.set_initial_time(dt0)
    # simplify by setting Greenwich phase to 0 and undoing nodal corrections
    tide.phi[:] = 0.
    tide.f[:] = 1
    tide.u[:] = 0.
    tnci = uptide.TPXOTidalInterpolator(tide, dummy_tpxo_grid_file, dummy_tpxo_elev_file)

    for t in [0., 1000., 86400.]:
        tnci.set_time(t)
        np.testing.assert_allclose(tnci.get_val([0., 0.]), np.cos(tide.omega[0]*t))


def test_dummy_tpxo_multiple(dummy_tpxo_grid_file, dummy_tpxo_elev_file):
    selection = [0, 1, 4]
    tide = uptide.Tides([constituents[i] for i in selection])
    dt0 = datetime.datetime(2003, 3, 28, 0, 0, 0)
    tide.set_initial_time(dt0)
    # simplify by setting Greenwich phase to 0 and undoing nodal corrections
    tide.phi[:] = 0.
    tide.f[:] = 1
    tide.u[:] = 0.
    tnci = uptide.TPXOTidalInterpolator(tide, dummy_tpxo_grid_file, dummy_tpxo_elev_file)

    for t in [0., 1000., 86400.]:
        tnci.set_time(t)
        expected = sum([(i+1.0)*np.cos(omega*t) for i, omega in zip(selection, tide.omega)])
        np.testing.assert_allclose(tnci.get_val([0., 0.]), expected)


@pytest.fixture
def dummy_tpxo_atlas_files(tmp_path):
    """Dummy files in TPXO atlas layout: 1D coordinates, a grid file with bathymetry hz
    (zero on land) but no mask, and a separate elevation file per constituent with
    integer hRe and hIm in millimeters."""
    nx, ny = 4, 3
    lon = np.arange(nx) * 90. + 45.  # fully periodic: 45, 135, 225, 315
    lat = np.array([-45., 0., 45.])
    ds = netCDF4.Dataset(tmp_path / 'grid_atlas.nc', 'w', format='NETCDF3_CLASSIC')
    ds.createDimension('nx', nx)
    ds.createDimension('ny', ny)
    for c in 'zuv':
        ds.createVariable(f'lon_{c}', np.float64, ('nx',))[:] = lon
        ds.createVariable(f'lat_{c}', np.float64, ('ny',))[:] = lat
    hz = np.full((nx, ny), 100., dtype=np.float32)
    hz[2, :] = 0.  # land at lon=225
    for c in 'zuv':
        ds.createVariable(f'h{c}', np.float32, ('nx', 'ny'))[:] = hz
    ds.close()

    for i, constituent in enumerate(constituents):
        ds = netCDF4.Dataset(tmp_path / f'h_{constituent.lower()}_atlas.nc', 'w', format='NETCDF3_CLASSIC')
        ds.createDimension('nx', nx)
        ds.createDimension('ny', ny)
        nct = ds.createDimension('nct', 4)
        ds.createVariable('con', 'c', ('nct',))[:] = netCDF4.stringtochar(np.array([constituent.lower().ljust(nct.size)], 'S4'))[0]
        ds.createVariable('lon_z', np.float64, ('nx',))[:] = lon
        ds.createVariable('lat_z', np.float64, ('ny',))[:] = lat
        # give M2 amplitude of 1m, S2 amplitude of 2m, etc. - except on land
        val = np.full((nx, ny), (i+1)*1000, dtype=np.int32)
        val[2, :] = 0
        hRe = ds.createVariable('hRe', np.int32, ('nx', 'ny'))
        hRe.units = 'millimeter'
        hRe[:] = val
        hIm = ds.createVariable('hIm', np.int32, ('nx', 'ny'))
        hIm.units = 'millimeter'
        hIm[:] = 0
        ds.close()

    for i, constituent in enumerate(constituents):
        ds = netCDF4.Dataset(tmp_path / f'u_{constituent.lower()}_atlas.nc', 'w', format='NETCDF3_CLASSIC')
        ds.createDimension('nx', nx)
        ds.createDimension('ny', ny)
        nct = ds.createDimension('nct', 4)
        ds.createVariable('con', 'c', ('nct',))[:] = netCDF4.stringtochar(np.array([constituent.lower().ljust(nct.size)], 'S4'))[0]
        for c in 'uv':
            ds.createVariable(f'lon_{c}', np.float64, ('nx',))[:] = lon
            ds.createVariable(f'lat_{c}', np.float64, ('ny',))[:] = lat
        # transport amplitudes in cm^2/s: (i+1) m^2/s for u, and -(i+1) m^2/s for v,
        # so with a depth of 100m, velocity (i+1)/100 m/s for u and -(i+1)/100 m/s for v
        for c, sign in (('u', 1), ('v', -1)):
            val = np.full((nx, ny), sign*(i+1)*10000, dtype=np.int32)
            val[2, :] = 0
            re = ds.createVariable(f'{c}Re', np.int32, ('nx', 'ny'))
            re.units = 'centimeter^2/sec'
            re[:] = val
            im = ds.createVariable(f'{c}Im', np.int32, ('nx', 'ny'))
            im.units = 'centimeter^2/sec'
            im[:] = 0
        ds.close()

    yield str(tmp_path / 'grid_atlas.nc'), str(tmp_path / 'h_{con}_atlas.nc')


def _simple_tide(selection):
    tide = uptide.Tides([constituents[i] for i in selection])
    tide.set_initial_time(datetime.datetime(2003, 3, 28, 0, 0, 0))
    # simplify by setting Greenwich phase to 0 and undoing nodal corrections
    tide.phi[:] = 0.
    tide.f[:] = 1
    tide.u[:] = 0.
    return tide


@pytest.mark.parametrize('ranges', [None, ((-50., 50.), (-10., 10.))])
def test_dummy_tpxo_atlas(dummy_tpxo_atlas_files, ranges):
    grid_file, data_file = dummy_tpxo_atlas_files
    selection = [0, 1, 4]
    tide = _simple_tide(selection)
    tnci = uptide.TPXOTidalInterpolator(tide, grid_file, data_file, ranges=ranges)

    for t in [0., 1000., 86400.]:
        tnci.set_time(t)
        expected = sum([(i+1.0)*np.cos(omega*t) for i, omega in zip(selection, tide.omega)])
        # interpolation across the periodic boundary (between lon=315 and lon=45)
        np.testing.assert_allclose(tnci.get_val([0., 0.]), expected)
        np.testing.assert_allclose(tnci.get_val([-10., 20.]), expected)
        # next to land (lon=225), the land value (0) should be ignored
        np.testing.assert_allclose(tnci.get_val([150., 0.]), expected)


def test_dummy_tpxo_atlas_land(dummy_tpxo_atlas_files):
    grid_file, data_file = dummy_tpxo_atlas_files
    tide = _simple_tide([0])
    tnci = uptide.TPXOTidalInterpolator(tide, grid_file, data_file)
    tnci.set_time(0.)
    with pytest.raises(uptide.netcdf_reader.CoordinateError):
        tnci.get_val([225., 0.])


def test_dummy_tpxo_atlas_missing_constituent(dummy_tpxo_atlas_files):
    grid_file, data_file = dummy_tpxo_atlas_files
    tide = uptide.Tides(['M2', 'M4'])
    with pytest.raises(FileNotFoundError, match='M4'):
        uptide.TPXOTidalInterpolator(tide, grid_file, data_file)


def test_dummy_tpxo_atlas_wrong_constituent(dummy_tpxo_atlas_files, tmp_path):
    grid_file, data_file = dummy_tpxo_atlas_files
    os.rename(data_file.format(con='s2'), data_file.format(con='m4'))
    tide = uptide.Tides(['M4'])
    with pytest.raises(ValueError, match='contains constituent s2'):
        uptide.TPXOTidalInterpolator(tide, grid_file, data_file)


@pytest.fixture
def dummy_tpxo_uv_file(tmp_path):
    ds = netCDF4.Dataset(tmp_path / 'u_tpxo9.nc', 'w', format='NETCDF3_CLASSIC')
    ds.createDimension('nx', 2)
    ds.createDimension('ny', 2)
    ds.createDimension('nc', len(constituents))
    nct = ds.createDimension('nct', 4)
    con = ds.createVariable('con', 'c', ('nc', 'nct'))
    con[:] = [x.ljust(nct.size) for x in constituents]
    amp = (np.arange(len(constituents))[:, np.newaxis, np.newaxis] + 1.)*np.ones((1, 2, 2))
    for c in 'uv':
        # velocity amplitude in cm/s: M2 1 cm/s, S2 2 cm/s, etc.
        a = ds.createVariable(f'{c}a', np.float64, ('nc', 'nx', 'ny'))
        a.units = 'cm/s'
        a[:] = amp
        p = ds.createVariable(f'{c}p', np.float64, ('nc', 'nx', 'ny'))
        p.units = 'degree GMT'
        p[:] = 0.
        # transport amplitude in m^2/s: M2 10 m^2/s, S2 20 m^2/s, etc.
        A = ds.createVariable(f'{c.upper()}a', np.float64, ('nc', 'nx', 'ny'))
        A.units = 'meter^2/s'
        A[:] = 10*amp
        re = ds.createVariable(f'{c.upper()}Re', np.float64, ('nc', 'nx', 'ny'))
        re.units = 'meter^2/s'
        re[:] = 10*amp
        im = ds.createVariable(f'{c.upper()}Im', np.float64, ('nc', 'nx', 'ny'))
        im.units = 'meter^2/s'
        im[:] = 0.
    filepath = ds.filepath()
    ds.close()
    yield filepath
    os.remove(filepath)


@pytest.mark.parametrize('component', ['u', 'v'])
def test_dummy_tpxo_velocity(dummy_tpxo_grid_file, dummy_tpxo_uv_file, component):
    selection = [0, 1, 4]
    tide = _simple_tide(selection)
    # velocities are returned in m/s
    tnci = TPXOncTidalComponentInterpolator(tide, dummy_tpxo_grid_file, dummy_tpxo_uv_file,
                                            component, component)
    # transports in m^2/s
    tnciT = TPXOncTidalComponentInterpolator(tide, dummy_tpxo_grid_file, dummy_tpxo_uv_file,
                                             component, component.upper())
    for t in [0., 1000., 86400.]:
        tnci.set_time(t)
        tnciT.set_time(t)
        expected = sum([(i+1.0)*np.cos(omega*t) for i, omega in zip(selection, tide.omega)])
        np.testing.assert_allclose(tnci.get_val([0., 0.]), 0.01*expected)
        np.testing.assert_allclose(tnciT.get_val([0., 0.]), 10*expected)


@pytest.mark.parametrize('ranges', [None, ((-50., 50.), (-10., 10.))])
def test_dummy_tpxo_atlas_velocity(dummy_tpxo_atlas_files, ranges):
    grid_file, elev_file = dummy_tpxo_atlas_files
    uv_file = elev_file.replace('h_{con}', 'u_{con}')
    selection = [0, 1, 4]
    tide = _simple_tide(selection)
    for component, sign in (('u', 1), ('v', -1)):
        # velocities computed from transport/depth, in m/s
        tnci = TPXOncTidalComponentInterpolator(tide, grid_file, uv_file, component, component, ranges=ranges)
        # transports in m^2/s
        tnciT = TPXOncTidalComponentInterpolator(tide, grid_file, uv_file, component, component.upper(),
                                                 ranges=ranges)
        for t in [0., 1000., 86400.]:
            tnci.set_time(t)
            tnciT.set_time(t)
            expected = sign*sum([(i+1.0)*np.cos(omega*t) for i, omega in zip(selection, tide.omega)])
            for x in ([0., 0.], [150., 0.]):
                np.testing.assert_allclose(tnci.get_val(x), expected/100.)
                np.testing.assert_allclose(tnciT.get_val(x), expected)


def test_dummy_tpxo_velocity_missing(dummy_tpxo_atlas_files):
    grid_file, elev_file = dummy_tpxo_atlas_files
    tide = _simple_tide([0])
    # the elevation files do not contain velocities or transports
    with pytest.raises(ValueError, match='No velocity data'):
        TPXOncTidalComponentInterpolator(tide, grid_file, elev_file, 'u', 'u')
