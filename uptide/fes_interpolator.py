from six import with_metaclass
import abc
import datetime
import tempfile
import os
from contextlib import contextmanager
from .tides import Tides
import numpy


ALL_FES2014_TIDAL_CONSTITUENTS = ["2N2", "EPS2", "J1", "K1", "K2", "L2", "M2", "M3", "M4", "M6", "M8",
                                  "MF", "MKS2", "MM", "MN4", "MS4", "MSF", "MSQM", "MTM", "MU2", "N2", "N4", "NU2",
                                  "O1", "P1", "Q1", "R2", "S1", "S2", "S4", "SA", "SSA", "T2"]


class TidalInterpolator(with_metaclass(abc.ABCMeta)):
    """Abstract base class for tidal interpolators."""

    def set_initial_time(self, datetime0):
        """Set datetime corresponding to t=0

        If datetime without timezone info is supplied (tzinfo), datetime is assumed to be UTC."""
        if datetime0.tzinfo:
            # only pull in this dependency when tzinfo is supplied
            import pytz
            self.datetime0 = pytz.utc.localize(datetime0)
        else:
            # with naive datetimes without tzinfo, the assumption is everything is in UTC
            self.datetime0 = datetime0

    @abc.abstractmethod
    def set_time(self, t):
        pass

    @abc.abstractmethod
    def get_val(self, x, **kwargs):
        pass


fes_ini_template = """TIDE_{constituent}_FILE         = {fes_data_path}/{lower_case_constituent}.nc
TIDE_{constituent}_LATITUDE     = lat
TIDE_{constituent}_LONGITUDE    = lon
TIDE_{constituent}_AMPLITUDE    = amplitude
TIDE_{constituent}_PHASE        = phase

"""


@contextmanager
def temporary_fes_ini_file(tide, fes_data_path):
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as f:
        for constituent in tide.constituents:
            f.write(fes_ini_template.format(
                constituent=constituent, fes_data_path=fes_data_path,
                lower_case_constituent=constituent.lower()
            ))
    f.close()
    yield f
    file_name = f.name

    os.remove(file_name)


class FES2014TidalInterpolator(TidalInterpolator):
    """Tidal interpolator based on FES2014 global solution.

    For any given time and lat, lon, interpolates amplitude and phases of
    the harmonic constituents of the FES2014 global tide solution and reconstructs
    the tidal elevation in that point.

    Data can be downloaded from ftp://ftp.aviso.altimetry.fr/auxiliary/tide_model/fes2014_elevations_and_load/fes2014b_elevations/
    Download either ocean_tide.tar.xz or ocean_tide_extrapolated.tar.xz (latter extrapolates amplitudes and phases inland
    so that interpolation in coastal locations is less dependent on being exactly within the FES2014 grid) and extract
    the netcdf files. For ftp access yYou need prior registration as described here:
        https://www.aviso.altimetry.fr/en/data/products/auxiliary-products/global-tide-fes.html.

    To make use of this class you need to install the fes python package using:
        pip install git+https://bitbucket.org/cnes_aviso/fes.git#subdirectory=python.

    Constituents and initial (t=0) datetime are set via a Tides object. To use all constituents available
    in the FES2014 solution:

        import uptide
        import datetime

        tide = uptide.Tides(uptide.ALL_FES2014_TIDAL_CONSTITUENTS)
        tide.set_initial_time(datetime.datetime(1975, 10, 24, 0, 0))
        tnci = uptide.FES2014TidalInterpolator(tide, 'path_to_extracted_fes2014_solution/ocean_tide/')

    Alternatively, you can also use the .ini files as provided in the fes source code (https://bitbucket.org/cnes_aviso/fes),
    and set the initial time seperately:

        tnci = uptide.FES2014TidalInterpolator('path_to/ocean_tide.ini')
        tnci.set_initial_time(datetime.datetime(2005,3,1,0,0))

    The tidal elevation at time t (in seconds) is obtained via:

        tnci.set_time(self, t)
        eta = tnci.get_val(self, (lat, lon))

    Here -90<lat<90 and 0<lon<360. Finally, the long period (longer than a year) tidal constituents can be excluded via:

        tnci = uptide.FES2014TidalInterpolator(..., include_long_period=False)"""
    def __init__(self, tide_or_fes_ini_file,
                 fes_data_path=None, include_long_period=True):
        # only import here to avoid hard dependency on fes
        # there are two versions of this, the old (pre 2.9.1) is imported as fes,
        # but from 2.9.1 we need to `import pyfes`
        try:
            import fes
        except ImportError:
            try:
                import pyfes as fes
            except ImportError:
                raise ImportError("Failed to import fes. See https://github.com/stephankramer/uptide for installation instructions.")

        if isinstance(tide_or_fes_ini_file, Tides):
            self.set_initial_time(tide_or_fes_ini_file.datetime0)
            with temporary_fes_ini_file(tide_or_fes_ini_file, fes_data_path) as f:
                self.fh = fes.Handler("ocean", "io", f.name)
        else:
            assert fes_data_path is None, "Do not provide fes_data_path if fes_ini_file is specified"
            self.fh = fes.Handler("ocean", "io", tide_or_fes_ini_file)

        self.include_long_period = include_long_period

    def set_time(self, t):
        """Set time (in seconds) at which to reconstruct tide

        Time is with respect to datetime set via set_initial_time() method on the FES2014TidalInterpolator itself
        or the initial Tides object provided."""
        self.current_datetime = self.datetime0 + datetime.timedelta(seconds=t)

    def get_val(self, x):
        """Evaluate tide in location x=(lat, lon)

        Here -90<lat<90 and 0<lon<360."""
        try:
            # old API:
            st, lt = self.fh.scalar(x[0], x[1], self.current_datetime)
        except AttributeError:
            # new API
            st, lt, fes_min = self.fh.calculate(numpy.array([x[1]]), numpy.array([x[0]]), numpy.array([self.current_datetime]))
        # FES2014 is in cm, others are all in m
        if self.include_long_period:
            return (st+lt) * 0.01
        else:
            return st * 0.01
