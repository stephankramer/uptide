About uptide
==============
uptide is a python package for tidal calculations.

Prerequisites
---------------
* python 2.6 or 2.7 (not tested for python 3)
* numpy
* optional: scipy with netCDF support (scipy.io.netcdf)

Functionality
---------------
* Given the phase and amplitudes of the harmonic constituents (M2, S2, etc.) reconstruct the tidal signal at an arbitrary date and time (including nodal corrections).
* Read the phases and amplitudes on a regular Cartesian grid from a netCDF file and interpolate the tidal signal at any arbitrary point in the domain, for any date and time. Phases and amplitudes can be read directly from netCDF files as provided by the [FES2004][1] and [OTPSnc][2] global and regional tidal data bases.
[1]:ftp://ftp.legos.obs-mip.fr/pub/soa/maree/tide_model/global_solution/fes2004/
[2]:http://volkov.oce.orst.edu/tides/otps.html
