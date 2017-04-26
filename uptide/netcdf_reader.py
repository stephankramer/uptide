# horrible kludge to import python netcdf class - there's three different implementations to choose from!
# luckily they adhere to the same API
try:
  # this seems to be the most mature and also handles netcdf4 files
  from netCDF4 import Dataset as NetCDFFile
except ImportError:
  try:
    # this one is older but is quite often installed
    from Scientific.IO.NetCDF import NetCDFFile
  except ImportError:
    # finally, try the one in scipy that I hear conflicting things about
    try:
      # this only works in python 2.7
      from scipy.io.netcdf import NetCDFFile
    except ImportError:
      # in python 2.6 it's called something else
      from scipy.io.netcdf import netcdf_file as NetCDFFile
import math
import numpy
import numpy.ma

# any error generated by the NetCDFInterpolator object:
class NetCDFInterpolatorError(Exception): pass

# error caused by coordinates being out of range or inside land mask:
class CoordinateError(NetCDFInterpolatorError):
  def __init__(self, message, x, i, j):
    self.message = message
    self.x = x
    self.ij = i, j
  def __str__(self):
    return "at x, y={} indexed at i, j={}; {}".format(self.x, self.ij, self.message)

class Interpolator(object):
  def __init__(self, origin, delta, val, mask=None):
    self.origin = origin
    self.delta = delta
    self.val = val
    self.mask = mask
    # cache points that need to be extrapolated
    self.extrapolation_points = {}

  def set_mask(self, mask):
    self.mask = mask
    # changing the mask invalidates the extrapolation cache
    self.extrapolation_points ={}

  def find_extrapolation_points(self, x, i, j):
    if x in self.extrapolation_points:
      return self.extrapolation_points[x]

    # This should only happen infrequently, so warn user (commented out because a user complained):
    # print "Need to extrapolate point coordinates ", x
    ijs = [(i-1, j+1), (i-1, j), (i, j-1), (i+1, j-1), (i+2, j), (i+2, j+1), (i+1, j+2), (i, j+2)] # Neighbouring points
    ijs += [(i-1, j-1), (i+2, j-1), (i+2, j+2), (i-1, j+2)] # Diagonal points

    extrap_points = []
    for a, b in ijs:
      try:
        if self.mask[a, b]:
          extrap_points.append((a, b))
      except IndexError:
        # if we go out of range, ignore this point
        pass

    if len(extrap_points) == 0:
      raise CoordinateError("Inside landmask - tried extrapolating but failed", x, i, j)

    self.extrapolation_points[x] = extrap_points

    return extrap_points


  def get_val(self, x, allow_extrapolation=False):
    xhat = (x[0]-self.origin[0])/self.delta[0]
    yhat = (x[1]-self.origin[1])/self.delta[1]
    i = int(math.floor(xhat))
    j = int(math.floor(yhat))
    # this is not catched as an IndexError below, because of wrapping of negative indices
    if i<0 or j<0:
      raise CoordinateError("Coordinate out of range", x, i, j)
    alpha = xhat % 1.0
    beta = yhat % 1.0
    try:
      if self.mask is not None:

        # case with a land mask

        w00 = (1.0-alpha)*(1.0-beta)*self.mask[i,j]
        w10 = alpha*(1.0-beta)*self.mask[i+1,j]
        w01 = (1.0-alpha)*beta*self.mask[i,j+1]
        w11 = alpha*beta*self.mask[i+1,j+1]
        if len(self.val.shape)==2:
          value = w00*self.val[i,j] + w10*self.val[i+1,j] + w01*self.val[i,j+1] + w11*self.val[i+1,j+1]
        elif len(self.val.shape)==3:
          value = w00*self.val[:,i,j] + w10*self.val[:,i+1,j] + w01*self.val[:,i,j+1] + w11*self.val[:,i+1,j+1]
        else:
          raise NetCDFInterpolatorError("Field to interpolate, should have 2 or 3 dimensions")
        sumw = w00+w10+w01+w11

        if sumw>0.0:
          value = value/sumw
        elif allow_extrapolation:
          extrap_points = self.find_extrapolation_points(x, i, j)
          if len(self.val.shape)==2:
            value = sum([self.val[a, b] for a, b in extrap_points])/len(extrap_points)
          elif len(self.val.shape)==3:
            value = sum([self.val[:, a, b] for a, b in extrap_points])/len(extrap_points)
          else:
            raise NetCDFInterpolatorError("Field to interpolate, should have 2 or 3 dimensions")
        else:
          raise CoordinateError("Probing point inside land mask", x, i, j)

      else:

        # case without a land mask

        if len(self.val.shape)==2:
          value = ((1.0-beta)*((1.0-alpha)*self.val[i,j]+alpha*self.val[i+1,j])+
                  beta*((1.0-alpha)*self.val[i,j+1]+alpha*self.val[i+1,j+1]))
        elif len(self.val.shape)==3:
          value = ((1.0-beta)*((1.0-alpha)*self.val[:,i,j]+alpha*self.val[:,i+1,j])+
                  beta*((1.0-alpha)*self.val[:,i,j+1]+alpha*self.val[:,i+1,j+1]))
        else:
          raise NetCDFInterpolatorError("Field to interpolate, should have 2 or 3 dimensions")
    except IndexError:
      raise CoordinateError("Coordinate out of range", x, i, j)
    return value
    
# note that a NetCDFInterpolator is *not* object an Interpolator object
# the latter is considered immutable, whereas the NetCDFInterpolator may
# change with set_ranges() and set_field() and will create a new Interpolator sub-object
# each time
class NetCDFInterpolator(object):
  """Implements an object to interpolate values from a NetCDF-stored data set.

  The NetCDF file should contain two coordinate fields, e.g. latitude and longitude. Each of those two coordinates
  is assumed to be aligned with one dimension of the logical 2D grid and should be equi-distant. 
  To open the NetCDFInterpolator object:

    nci = NetCDFInterpolator('foo.nc', ('nx', 'ny'), ('longitude', latitude'))

  Here 'nx' and 'ny' refer to the names of the dimensions of the logical 2D grid and 'longitude' and 'latitude' to the 
  names of the coordinate fields. The order of these should match (e.g. here, the 'nx' dimension should be the 
  dimension in which 'longitude' increases, and 'ny' the dimension in which 'latitude' increases) and determines the order
  of coordinates in all other arguments (ranges and get_val).
  The names  can of dimension and coordinate fields can be obtained by using the ncdump program of the standard NetCDF utils:
    
    $ ncdump -h foo.nc
    netcdf foo {
      dimensions:
        nx = 20 ;
        ny = 10 ;
        variables:
        double z(nx, ny) ;
        double mask(nx, ny) ;
        double longitude(nx) ;
        double latitude(ny) ;
    }

  The coordinate fields may be stored as 1d or 2d fields (although only two values will actually be read to determine the 
  origin and step size). The order of the dimensions and coordinate fields specified in the call does not have to match that
  of the netCDF file, i.e. we could have opened the same file with:

    nci_transpose = NetCDFInterpolator('foo.nc', ('ny', 'nx'), ('latitude', longitude'))

  The only difference would be the order in which the coordinates to get_val and set_ranges are specified.
  To indicate the field to be interpolated:

    nci.set_field('z')

  To interpolate this field in any arbitrary point:

    nci.get_val((-3.0, 58.5))

  If many interpolations are done 
  within a sub-domain of the area covered by the NetCDF, it may be much more efficient to indicate the range of coordinates
  with:

     nci.set_ranges(((-4.0,-2.0),(58.0,59.0)))

  This will load all values within the indicated range (here -4.0<longitude<-2.0 and 58.0<latitude<59.0) in memory.
  A land-mask can be provided to avoid interpolating from undefined land-values. The mask field should be 0.0 in land points
  and 1.0 at sea.

     nci.set_mask('mask')

  Alternatively, a mask can be defined from a fill value that has been used to indicate undefined land-points. The field name
  (not necessarily the same as the interpolated field) and fill value should be provided:

     nci.set_mask_from_fill_value('z', -9999.)

  It is allowed to switch between different fields using multiple calls of set_field(). The mask and ranges will be retained. It is
  however not allowed to call set_mask() or set_ranges() more than once. Finally, 
  for the case where the coordinate fields (and optionally
  the mask field) is stored in a different file than the one containing the field values to be interpolated, the following syntax
  is provided:

    nci1 = NetCDFInterpolator('grid.nc', ('nx', 'ny'), ('longitude', latitude'))
    nci1.set_mask('mask')
    nci2 = NetCDFInterpolator('values.nc', nci1)
    nci2.set_field('temperature')
    nci2.get_val(...)

  Here, the coordinate information of nci, including the mask and ranges if set, are copied and used in nci2.

  """
  def __init__(self, filename, *args, **kwargs):
    self.nc = NetCDFFile(filename, 'r')

    if len(args)==1:

      # we copy the grid information of another netcdf interpolator

      nci = args[0]
      self.dimensions = nci.dimensions
      self.shape = nci.shape
      self.origin = nci.origin
      self.delta = nci.delta
      self.iranges = nci.iranges
      self.mask = nci.mask
      if nci.mask is not None:
        self.dim_order = nci.dim_order

    elif len(args)==2:

      dimensions = args[0]
      coordinate_fields = args[1]

      self.dimensions = dimensions
      self.shape = []
      self.origin = []
      self.delta = []

      for dimension,field_name in zip(dimensions, coordinate_fields):
        N = self.nc.dimensions[dimension]
        if not isinstance(N, int):
          # let's guess it's a netCDF4.Dimension, so we should ask for its len (yuck)
          N = len(N)
        self.shape.append(N)
        val = self.nc.variables[field_name]
        if len(val.shape)==1:
          self.origin.append(val[0])
          self.delta.append((val[-1]-self.origin[-1])/(N-1))
        elif len(val.shape)==2:
          self.origin.append(val[0,0])
          self.delta.append((val[-1,-1]-self.origin[-1])/(N-1))
        else:
          raise NetCDFInterpolatorError("Unrecognized shape of coordinate field")

      self.iranges = None
      self.mask = None

    self.interpolator = None

    if "ranges" in kwargs:
      ranges = kwargs("ranges")
      self.set_ranges(ranges)

  def set_ranges(self, ranges):
    """Set the range of the coordinates. All the values of points located within this range are read from file at once.
    This may be more efficient if many interpolations are done within this domain."""
    if self.iranges is not None:
      # this could probably be fixed, but requires thought and testing:
      raise NetCDFInterpolatorError("set_ranges() should only be called once!")

    self.iranges = []
    origin_new = []
    shape_new = []
    for xlimits,xmin,xshape,deltax in zip(ranges,self.origin,self.shape,self.delta):
      # compute the index range imin:imax
      # for the min, take one off (i.e. add an extra row column) to avoid rounding issues:
      imin = max( int((xlimits[0]-xmin)/deltax)-1, 0 )
      # for the max, we add 3 because:
      # 1) we're rounding off first
      # 2) add one extra for rounding off issues
      # 3) python imin:imax range means up to and *excluding* imax
      # Example: xmin=0.0,deltax=1.0,xlimits[1]=3.7 -> imax=6 
      # (which is one too many, but nearing 3.999 we may get into a situation where we have to interpolate between 4 and 5)
      imax = min( int((xlimits[1]-xmin)/deltax)+3, xshape)
      if imin>=imax:
        raise NetCDFInterpolatorError("Provided ranges outside netCDF range")
      self.iranges.append((imin,imax))
      origin_new.append(xmin+imin*deltax)
      shape_new.append(imax-imin)

    self.origin = origin_new
    self.shape = shape_new

    if self.mask is not None:
      ir = [self.iranges[d] for d in self.dim_order]
      self.mask = self.mask[ir[0][0]:ir[0][1],ir[1][0]:ir[1][1]]

    if self.interpolator is not None:
      ir = [self.iranges[d] for d in self.dim_order]
      if len(self.val.shape)==2:
        self.val = self.val[ir[0][0]:ir[0][1],ir[1][0]:ir[1][1]]
      elif len(self.val.shape)==3:
        self.val = self.val[:,ir[0][0]:ir[0][1],ir[1][0]:ir[1][1]]
      origin = [self.origin[d] for d in self.dim_order]
      delta = [self.delta[d] for d in self.dim_order]
      self.interpolator = Interpolator(origin, delta, self.val, self.mask)

  def set_mask(self, field_name):
    """Sets a land mask from a mask field. This field should have a value of 0.0 for land points and 1.0 for the sea"""
    mask = self.nc.variables[field_name]
    if hasattr(mask, 'set_auto_maskandscale'):
      # netCDF4 automagically turns a variable into a masked array if it has a _FillValue attribute
      # might be handy, but is not compatible with other netcdf implementations
      mask.set_auto_maskandscale(False)

    # work out the dimension, in particular its order
    if list(mask.dimensions) == list(self.dimensions):
      dim_order = [0,1]
    elif list(mask.dimensions) == list(self.dimensions)[::-1]:
      dim_order = [1,0]
    else:
      raise NetCDFInterpolatorError("Dimensions of mask field not the same as specified in __init__")

    if self.iranges is not None:
      ir = [self.iranges[d] for d in dim_order]
      mask = mask[ir[0][0]:ir[0][1],ir[1][0]:ir[1][1]]

    self._set_mask_and_dim_order(mask, dim_order)

  def _set_mask_and_dim_order(self, mask, mask_dim_order):
    """If no interpolator is defined yet, sets dim_order to mask_dim_order otherwise
    ensures that the dim_order of the interpolator is the same as that of mask_dim_order, if
    necessary by transposing the mask."""

    if self.interpolator is None:
      self.mask = mask
      # dim_order of mask is always the same as that of the field (mask will be transposed if the next field is in different order)
      self.dim_order = mask_dim_order
    else:
      # dim_order is already set by interpolator, check that it agrees with mask_dim_order
      if mask_dim_order==self.dim_order:
        self.mask = mask
      else:
        if hasattr(mask, 'T'):
          # if we've made a copy into a numpy array already, we can do a cheap in-place transpose
          self.mask = mask.T
        else:
          # requires a copy (and complete read from disk)
          self.mask = numpy.transpose(mask[:,:])
      self.interpolator.set_mask(self.mask)

  def set_mask_from_fill_value(self, field_name, fill_value):
    """Sets a land mask, where all points for which the supplied field equals the supplied fill value. The supplied field_name
    does not have to be the same as the field that is interpolated from, set with set_field()."""
    val = self.nc.variables[field_name]
    # work out the dimension, in particular its order
    if list(val.dimensions)[-2:] == list(self.dimensions):
      dim_order = [0,1]
    elif list(val.dimensions)[-2:] == list(self.dimensions)[::-1]:
      dim_order = [1,0]
    else:
      raise NetCDFInterpolatorError("Dimensions of mask field not the same as specified in __init__")

    if self.iranges is not None:
      ir = [self.iranges[d] for d in dim_order]
      if len(val.shape)==2:
        val = val[ir[0][0]:ir[0][1],ir[1][0]:ir[1][1]]
      elif len(val.shape)==3:
        # multiple values per gridpoint, just take the first one
        val = val[0,ir[0][0]:ir[0][1],ir[1][0]:ir[1][1]]
      else:
        raise NetCDFInterpolatorError("Field to extract mask from, should have 2 or 3 dimensions")
    elif len(val.shape)==2:
      val = val[:]
    elif len(val.shape)==3:
      # multiple values per gridpoint, just take the first one
      val = val[0,:,:]
    else:
      raise NetCDFInterpolatorError("Field to extract mask from, should have 2 or 3 dimensions")

    mask = numpy.logical_not(numpy.isclose(val, fill_value))
    self._set_mask_and_dim_order(mask, dim_order)

  def set_field(self, field_name):
    """Set the name of the field to be interpolated."""
    self.field_name = field_name
    self.val = self.nc.variables[field_name]

    # work out the dimension, in particular its order
    if list(self.val.dimensions)[-2:] == list(self.dimensions):
      dim_order = [0,1]
    elif list(self.val.dimensions)[-2:] == list(self.dimensions)[::-1]:
      dim_order = [1,0]
    else:
      raise NetCDFInterpolatorError("Dimensions of field not the same as specified in __init__")

    if self.iranges is not None:
      ir = [self.iranges[d] for d in dim_order]
      if len(self.val.shape)==2:
        self.val = self.val[ir[0][0]:ir[0][1],ir[1][0]:ir[1][1]]
      elif len(self.val.shape)==3:
        self.val = self.val[:,ir[0][0]:ir[0][1],ir[1][0]:ir[1][1]]
      else:
        raise NetCDFInterpolatorError("Field to interpolate, should have 2 or 3 dimensions")

    if self.mask is not None:
      if not self.dim_order==dim_order:
        if hasattr(self.mask, 'T'):
          # if we've made a copy into a numpy array already, we can do a cheap in-place transpose
          self.mask = self.mask.T
        else:
          # requires a copy (and complete read from disk)
          self.mask = numpy.transpose(self.mask[:,:])
    self.dim_order = dim_order

    origin = [self.origin[d] for d in self.dim_order]
    delta = [self.delta[d] for d in self.dim_order]
    self.interpolator = Interpolator(origin, delta, self.val, self.mask)

  def get_val(self, x, allow_extrapolation=False):
    """Interpolate the field chosen with set_field(). The order of the coordinates should correspond with the storage order in the file."""
    if not hasattr(self, "interpolator"):
      raise NetCDFInterpolatorError("Should call set_field() before calling get_val()!")
    if self.dim_order[0]==0:
      return self.interpolator.get_val(x, allow_extrapolation)
    else:
      # swap dimensions
      return self.interpolator.get_val((x[1],x[0]), allow_extrapolation)
