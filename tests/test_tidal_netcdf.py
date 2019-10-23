import pytest
import os
pytestmark = pytest.mark.skipif(
    not any(x in os.environ for x in ('AMCG_TIDAL_FILE', 'OTPS_GRID_FILE', 'TPXO_GRID_FILE',
                                      'FES2004_FILE', 'FES2014_INI_FILE', 'FES2014_DATA_PATH')),
    reason="Only run when tidal data file are given with environment variables"
)

import uptide.tidal_netcdf
import uptide
import datetime
import numpy
from numpy.testing import assert_equal

try:
    import matplotlib.pyplot as plt
except ImportError:
    plt = None


def extract_series(tnci, ll, trange):
    etas = []
    for t in trange:
        tnci.set_time(t)
        eta = tnci.get_val(ll)
        etas.append(eta)
    return etas


def test():
    constituents = ('M2', 'S2', 'N2', 'K2', 'K1', 'O1', 'P1', 'Q1')
    lat = 58.708
    lon = -3.287
    lon = -1.
    dt0 = datetime.datetime(2003, 3, 28, 0, 0, 0)
    trange = numpy.arange(0., 24.*32., 0.25)*3600.

    tide = uptide.Tides(constituents)
    tide.set_initial_time(dt0)

    series = {}

    if 'AMCG_TIDAL_FILE' in os.environ:
        tnci = uptide.tidal_netcdf.AMCGTidalInterpolator(tide,
                                                         os.environ['AMCG_TIDAL_FILE'],
                                                         ranges=((58.0, 61.0), (-4.0, 0.0)))
        series['AMCG'] = extract_series(tnci, (lat, lon), trange)

    if 'OTPS_GRID_FILE' in os.environ:
        tnci = uptide.tidal_netcdf.TPXOTidalInterpolator(tide,
                                                         os.environ['OTPS_GRID_FILE'],
                                                         os.environ['OTPS_DATA_FILE'],
                                                         ranges=((-4.0, 0.0), (58.0, 61.0)))
        series['OTPS'] = extract_series(tnci, (lon, lat), trange)
    if 'TPXO_GRID_FILE' in os.environ:
        tnci = uptide.tidal_netcdf.TPXOTidalInterpolator(tide,
                                                         os.environ['TPXO_GRID_FILE'],
                                                         os.environ['TPXO_DATA_FILE'],
                                                         ranges=((356., 360.0), (58.0, 61.0)))
        series['TXPO'] = extract_series(tnci, (360. + lon, lat), trange)

    if 'FES2004_FILE' in os.environ:
        tnci = uptide.tidal_netcdf.FESTidalInterpolator(tide,
                                                        os.environ['FES2004_FILE'],
                                                        ranges=((58.0, 61.0), (356., 360.0)))
        series['FES2004'] = extract_series(tnci, (lat, 360. + lon), trange)

    if 'FES2014_INI_FILE' in os.environ:
        tnci = uptide.FES2014TidalInterpolator(os.environ['FES2014_INI_FILE'])
        tnci.set_initial_time(dt0)
        series['FES2014'] = extract_series(tnci, (lat, 360. + lon), trange)

    if 'FES2014_INI_FILE' in os.environ:
        tnci = uptide.FES2014TidalInterpolator(os.environ['FES2014_INI_FILE'], include_long_period=False)
        tnci.set_initial_time(dt0)
        series['FES2014 short'] = extract_series(tnci, (lat, 360. + lon), trange)

    if 'FES2014_DATA_PATH' in os.environ:
        # check that the FES2014 tide is the same using a Tides object with all the
        # constituents that are in the standard ocean_tide.ini file
        tide2 = uptide.Tides(uptide.ALL_FES2014_TIDAL_CONSTITUENTS)
        tide2.set_initial_time(dt0)
        tnci = uptide.FES2014TidalInterpolator(tide2, os.environ['FES2014_DATA_PATH'])
        series['FES2014 all specified'] = extract_series(tnci, (lat, 360. + lon), trange)
        if 'FES2014' in series:
            assert_equal(series['FES2014'], series['FES2014 all specified'])

    if plt is None:
        raise Warning("Could not import matplotlib.pyplot")
    else:
        plt.figure()
        for label, etas in series.items():
            plt.plot(trange/uptide.tidal.day, etas, label=label)
        plt.legend()
        plt.show()
