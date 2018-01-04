About uptide
==============
uptide is a python package for tidal calculations. It computes tidal
free surface heights or velocities from the amplitudes and phases of the tidal
constituents. These amplitudes and phases can be read from global tidal
solutions such as [TPXO](http://volkov.oce.orst.edu/tides/) or [FES2014](https://www.aviso.altimetry.fr/en/data/products/auxiliary-products/global-tide-fes.html).
They can be read directly from the netCDF files provided by these sources. Some
limited functionality for tidal harmonic analysis is also available,

Prerequisites
---------------
* python 3 or 2.7 (deprecated)
* numpy
* to read from netCDF sources: python netCDF support. The
[netCDF4](https://github.com/Unidata/netcdf4-python) package is 
recommended. To install:
```
sudo CC=mpicc pip install netcdf4
```

or use the python-netcdf4 package on Ubuntu and Debian.
* for FES2014 support: the [FES package](https://bitbucket.org/cnes_aviso/fes). To install:
```
pip install git+https://bitbucket.org/cnes_aviso/fes.git#subdirectory=python
```

Functionality
---------------
* Given the phase and amplitudes of the harmonic constituents (M2, S2, etc.) reconstruct the tidal signal at an arbitrary date and time (including nodal corrections).
* Read the phases and amplitudes on a regular Cartesian grid from a netCDF file and interpolate the tidal signal at any arbitrary point in the domain, for any date and time. Phases and amplitudes can be read directly from netCDF files as provided by the [FES2004][1], [FES2014][2] and [TPXO][3] global and regional tidal data bases.

[1]: ftp://ftp.legos.obs-mip.fr/pub/soa/maree/tide_model/global_solution/fes2004/
[2]: https://www.aviso.altimetry.fr/en/data/products/auxiliary-products/global-tide-fes.html
[3]: http://volkov.oce.orst.edu/tides/otps.html
