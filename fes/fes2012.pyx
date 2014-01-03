cimport libfes

class FesException(Exception):
  pass

cdef class Fes(object):
  """Class to compute the tide using the fes library"""
  cdef libfes.FES _fes_handle
  def __cinit__(self, const char* const path, radial=False, in_memory=False):
    """Inititate fes class from the .ini file specified by the path argument:
      Fes(const char* const path, radial=False, in_memory=False)
    Optional arguments radial and in_memory maybe used respectively, 
    to compute the radial tide, and to read the entire global constituent 
    database in memory (as opposed to querying the NetCDF files for local data only)"""
    cdef libfes.fes_enum_access mode
    cdef libfes.fes_enum_tide_type tide
    if radial:
      tide = libfes.FES_RADIAL
    else:
      tide = libfes.FES_TIDE
    if in_memory:
      mode = libfes.FES_MEM
    else:
      mode = libfes.FES_IO

    if libfes.fes_new(&self._fes_handle, tide, mode, path):
      raise FesException(self.error())

  # doesn't do anythin currently, but provides docstring. Note that introspection of the arguments
  # doesn't work on c-extension methods, so we explicitly put it in the docstring.
  def __init__(self, const char* const path, radial=False, in_memory=False):
    """Inititate fes class from the .ini file specified by the path argument:
      Fes(const char* const path, radial=False, in_memory=False)
    Optional arguments radial and in_memory maybe used respectively, 
    to compute the radial tide, and to read the entire global constituent 
    database in memory (as opposed to querying the NetCDF files for local data only)"""
    pass

  def __dealloc__(self):
    """Deallocates the Fes object."""
    libfes.fes_delete(self._fes_handle)

  cdef core(self, const double lat, const double lon, const double time, 
      double *h, double *h_lp):
    """Directly calls fes_core to compute tide."""
    if libfes.fes_core(self._fes_handle, lat, lon, time, h, h_lp):
      raise FesException(self.error())

  def get_tide(self, lat, lon, time):
    """Computes the tide at the given location and time (sum of short and long period):
      h = <Fes-object>.get_tide(lat, lon, time)"""
    cdef double h, h_lp
    self.core(lat, lon, time, &h, &h_lp)
    return h+h_lp

  def get_short_and_long_period_tide(self, lat, lon, time):
    """Computes the tide at the given location and time, separated in a short-period and long-period
    contribution:
      h_sp, h_lp = <Fes-object>.get_short_and_long_period_tides(lat, lon, time)"""
    cdef double h, h_lp
    self.core(lat, lon, time, &h, &h_lp)
    return h, h_lp

  # return the error message
  def error(self):
    """If an error has occured this returns the error message produced by the fes library."""
    return libfes.fes_error(self._fes_handle)
  # return the error number
  def errno(self):
    """If an error has occured this returns the error number produced by the fes library."""
    return libfes.fes_errno(self._fes_handle)
