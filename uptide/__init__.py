"""
Python package for tidal computations. Main components:

    Tides - class of objects for tidal calculations. Provided with amplitudes and
        phases it computes the tidal signal

    harmonic_analysis   - function to compute harmonic analysis on a
        signal. Reconstructs amplitudes and phases of specified constituents.

    netcdf_reader.NetCDFInterpolator - class of objects to interpolate
        values stored in a NetCDF file.

    tidal_netcdf.TidalNetCDFInterpolator, tidal_netcdf.TPXOTidalInterpolator,
    tidal_netcdf.FESTidalInterpolator, tidal_netcdf.AMCGTidalInterpolator -
        class of objects to reconstruct tidal signal from tidal database of
        amplitudes and phases of tidal constituents stored in NetCDF format.
        Different formats, OTPSnc, FES and AMCG, are supported.

"""

from uptide.tides import Tides, select_constituents  # NOQA
from .analysis import harmonic_analysis  # NOQA
from .tidal_netcdf import OTPSncTidalInterpolator  # NOQA
from .tidal_netcdf import TPXOTidalInterpolator  # NOQA
from .fes_interpolator import FES2014TidalInterpolator, ALL_FES2014_TIDAL_CONSTITUENTS  # NOQA
from .ellipse import tidal_ellipse_parameters  # NOQA
