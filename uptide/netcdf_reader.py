from scipy.io.netcdf import NetCDFFile
import numpy

class Interpolator(object):
  def __init__(self, origin, delta, val, mask=None):
    self.origin = origin
    self.delta = delta
    self.val = val
    self.mask = mask

  def set_mask(self, mask):
    self.mask = mask

  def get_val(self, x):
    xhat = (x[0]-self.origin[0])/self.delta[0]
    yhat = (x[1]-self.origin[1])/self.delta[1]
    i = int(xhat)
    j = int(yhat)
    alpha = xhat % 1.0
    beta = yhat % 1.0
    try:
      if self.mask is not None:
        w00 = (1.0-alpha)*(1.0-beta)*self.mask[i,j]
        w10 = alpha*(1.0-beta)*self.mask[i+1,j]
        w01 = (1.0-alpha)*beta*self.mask[i,j+1]
        w11 = alpha*beta*self.mask[i+1,j+1]
        value = w00*self.val[...,i,j] + w10*self.val[...,i+1,j] + w01*self.val[...,i,j+1] + w11*self.val[...,i+1,j+1]
        sumw = w00+w10+w01+w11
        if sumw==0.0:
          raise Exception("Probing point inside land mask") 
        else:
          value = value/sumw
      else:
        value = ((1.0-beta)*((1.0-alpha)*self.val[...,i,j]+alpha*self.val[...,i+1,j])+
                  beta*((1.0-alpha)*self.val[...,i,j+1]+alpha*self.val[...,i+1,j+1]))
    except IndexError:
      print "Missing value at x,y,i,j = ", x[0], x[1], i, j
      raise
    return value
    
class NetCDFInterpolator(object):
  def __init__(self, filename, *args, **kwargs):
    self.nc = NetCDFFile(filename)

    if len(args)==1:

      # we copy the grid information of another netcdf interpolator

      nci = args[0]
      self.shape = nci.shape
      self.origin = nci.origin
      self.delta = nci.delta
      self.iranges = nci.iranges
      self.mask = nci.mask

    elif len(args)==2:

      dimensions = args[0]
      coordinate_fields = args[1]

      self.shape = []
      self.origin = []
      self.delta = []

      for dimension,field_name in zip(dimensions, coordinate_fields):
        self.shape.append(self.nc.dimensions[dimension])
        val = self.nc.variables[field_name]
        if len(val.shape)==1:
          self.origin.append(val[0])
          self.delta.append((val[-1]-self.origin[-1])/(self.shape[-1]-1))
        elif len(val.shape)==2:
          self.origin.append(val[0,0])
          self.delta.append((val[-1,-1]-self.origin[-1])/(self.shape[-1]-1))
        else:
          raise Exception("Unrecognized shape of coordinate field")

      self.iranges = None
      self.mask = None

    self.interpolator = None

    if "ranges" in kwargs:
      ranges = kwargs("ranges")
      self.set_ranges(ranges)

  def set_ranges(self, ranges):
    if self.iranges is not None:
      # this could probably be fixed, but requires thought and testing:
      raise Exception("set_ranges() should only be called once!")

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
        raise Exception("Provided ranges outside netCDF range")
      self.iranges.append((imin,imax))
      origin_new.append(xmin+imin*deltax)
      shape_new.append(imax-imin)

    self.origin = origin_new
    self.shape = shape_new

    if self.mask is not None:
      ir = self.iranges
      self.mask = self.mask[ir[0][0]:ir[0][1],ir[1][0]:ir[1][1]]
    if self.interpolator is not None:
      ir = self.iranges
      self.val = self.val[...,ir[0][0]:ir[0][1],ir[1][0]:ir[1][1]]
      self.interpolator = Interpolator(self.origin, self.delta, self.val, self.mask)

  def set_mask(self, field_name):
    if self.iranges is None:
      self.mask = self.nc.variables[field_name]
    else:
      ir = self.iranges
      self.mask = self.nc.variables[field_name][ir[0][0]:ir[0][1],ir[1][0]:ir[1][1]]
    if self.interpolator is not None:
      self.interpolator.set_mask(mask)

  def set_mask_from_fill_value(self, field_name, fill_value):
    if self.iranges is None:
      val = self.nc.variables[field_name]
    else:
      ir = self.iranges
      val = self.nc.variables[field_name][...,ir[0][0]:ir[0][1],ir[1][0]:ir[1][1]]
    if len(val.shape)==3:
      val = val[0,...]
    elif len(val.shape)==2:
      val = val[:]
    else:
      raise Exception("Field to extract mask from, should have 2 or 3 dimensions")
    self.mask = numpy.where(val==fill_value,0.,1.)

  def set_field(self, field_name):
    self.field_name = field_name
    if self.iranges is None:
      self.val = self.nc.variables[field_name]
    else:
      ir = self.iranges
      self.val = self.nc.variables[field_name][...,ir[0][0]:ir[0][1],ir[1][0]:ir[1][1]]
    self.interpolator = Interpolator(self.origin, self.delta, self.val, self.mask)


  def get_val(self, x):
    if not hasattr(self, "interpolator"):
      raise Exception("Should call set_field() before calling get_val()!")
    return self.interpolator.get_val(x)
