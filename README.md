About uptide
==============
uptide is a python package for tidal calculations. It computes the tidal
free surface height or velocities from the amplitudes and phases of the tidal
constituents. These amplitudes and phases can be read from global tidal
solutions such as [OSU/OTIS](http://volkov.oce.orst.edu/tides/) or [FES2012](http://www.aviso.oceanobs.com/en/data/products/auxiliary-products/global-tide-fes.html).
They can be read directly from the netCDF files provided by these sources. Some
limited functionality for tidal harmonic analysis is also available,

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
