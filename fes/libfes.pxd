cdef extern from "fes.h":
  ctypedef enum fes_enum_access:
    # Direct access (Grids are not loaded into memory).
    FES_IO = 0,
    # Memory access (Grids are loaded into memory).
    FES_MEM
  ctypedef enum fes_enum_tide_type:
    # Ocean tide
    FES_TIDE = 0,
    # Radial tide.
    FES_RADIAL
  ctypedef enum fes_enum_error:
    # No error reporting
    FES_SUCCESS,
    # Not enough memory
    FES_NO_MEMORY,
    # netCDF error
    FES_NETCDF_ERROR,
    # IO error
    FES_IO_ERROR,
    # Invalid configuration file
    FES_INI_ERROR,
    # No data available in grids for the location asked
    FES_NO_DATA
  ctypedef void* FES
  int fes_new(FES* handle, const fes_enum_tide_type tide,
        const fes_enum_access mode, const char* const path)
  void fes_delete(FES handle)
  int fes_core(FES handle, const double lat, const double lon, const double time,
        double* h, double* h_long_period)
  const char* fes_error(FES fes)
  fes_enum_error fes_errno(FES fes)
