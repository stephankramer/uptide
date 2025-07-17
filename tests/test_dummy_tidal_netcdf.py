import pytest
import uptide
import netCDF4
import os
import numpy as np
import datetime

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
