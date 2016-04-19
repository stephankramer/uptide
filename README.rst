About uptide
============

uptide is a python package for tidal calculations. It computes the tidal
free surface height or velocities from the amplitudes and phases of the
tidal constituents. These amplitudes and phases can be read from global
tidal solutions such as `OSU/OTIS <http://volkov.oce.orst.edu/tides/>`__
or
`FES2012 <http://www.aviso.oceanobs.com/en/data/products/auxiliary-products/global-tide-fes.html>`__.
They can be read directly from the netCDF files provided by these
sources. Some limited functionality for tidal harmonic analysis is also
available,

Prerequisites
-------------

-  python 2.6 or 2.7 (not tested for python 3)
-  numpy
-  to read from netCDF sources: python netCDF support. This can be
   either: `netCDF4 <http://code.google.com/p/netcdf4-python/>`__, or
   `Scientific.IO.NetCDF <http://dirac.cnrs-orleans.fr/plone/software/scientificpython/>`__,
   or `scipy.io.netcdf <http://www.scipy.org>`__. netCDF4 is the
   recommended package and is required for FES2012 that comes in netcdf4
   format. To install:

sudo CC=mpicc pip install netcdf4

For Ubuntu Precise, see this
`bug <http://code.google.com/p/netcdf4-python/issues/detail?id=194>`__.
So there, either use Scientific (sudo apt-get install
python-scientific), or install a newer version of netcdf4 (>=4.1.2).

Functionality
-------------

-  Given the phase and amplitudes of the harmonic constituents (M2, S2,
   etc.) reconstruct the tidal signal at an arbitrary date and time
   (including nodal corrections).
-  Read the phases and amplitudes on a regular Cartesian grid from a
   netCDF file and interpolate the tidal signal at any arbitrary point
   in the domain, for any date and time. Phases and amplitudes can be
   read directly from netCDF files as provided by the [FES2004][1],
   [FES2012][2] and [OTPSnc][3] global and regional tidal data bases.
   [1]:ftp://ftp.legos.obs-mip.fr/pub/soa/maree/tide\_model/global\_solution/fes2004/
   [2]:http://www.aviso.oceanobs.com/en/data/products/auxiliary-products/global-tide-fes.html
   [3]:http://volkov.oce.orst.edu/tides/otps.html
