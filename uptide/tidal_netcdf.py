import numpy
import netcdf_reader
import itertools

_deg2rad = numpy.pi/180.

class TidalNetCDFInterpolator(object):
  def __init__(self, tide, grid_file_name, dimensions, coordinate_fields,
      ranges=None, mask=None):
    self.tide = tide
    self.grid_file_name = grid_file_name
    self.nci = netcdf_reader.NetCDFInterpolator(grid_file_name, dimensions,
        coordinate_fields)

    if ranges is not None:
      self.set_ranges(ranges)
    if mask is not None:
      self.set_mask(mask)

  def set_ranges(self, ranges):
    self.nci.set_ranges(ranges)

  def set_mask(self, field_name):
    self.nci.set_mask(field_name)

  def set_mask_from_fill_value(self, field_name, fill_value):
    self.nci.set_mask_from_fill_value(field_name, fill_value)

  def load_amplitudes_and_phases(self, amplitude_file_name, amplitude_field_names,
      phase_file_name, phase_field_names):

    amp = numpy.array(self._collect_fields_val(amplitude_file_name, amplitude_field_names))
    phase = numpy.array(self._collect_fields_val(phase_file_name, phase_field_names))
    self.real_part = amp*numpy.cos(phase*_deg2rad)
    self.imag_part = -amp*numpy.sin(phase*_deg2rad)
    # NOTE: I did try several things with packing everything in a single nx x ny x nc x 2 array
    # (and making sure things are  contiguous in memory in the right order)
    # and contracting it with a nc x 2 array to compute the tides at all nx x ny points
    # but everything I tried was slower than this simple version

  def load_complex_components(self, real_file_name, real_field_names,
      imag_file_name, imag_field_names):
    self.real_part = numpy.array(self._collect_fields_val(real_file_name, real_field_names))
    self.imag_part = numpy.array(self._collect_fields_val(imag_file_name, imag_field_names))
    
  def _collect_fields_val(self, file_name, field_names):
    val = []
    if isinstance(file_name, basestring):
      file_names = itertools.repeat(file_name)
    else:
      file_names = file_name

    nci_filenm = self.grid_file_name
    nci = self.nci

    for filenm, fieldnm in zip(file_names, field_names):
      if filenm!=nci_filenm:
        # copies grid,mask and ranges information from self.nci (the "grid" netCDF file)
        nci = netcdf_reader.NetCDFInterpolator(filenm, self.nci)
        nci_filenm = filenm
      nci.set_field(fieldnm)
      val.append(nci.val[:])
    return val

  def load_amplitudes_and_phases_block(self,
      amplitude_file_name, amplitude_field_name, amplitude_field_components, 
      phase_file_name, phase_field_name, phase_field_components):
    
    amp = numpy.array(self._collect_fields_block(amplitude_file_name, amplitude_field_name, amplitude_field_components))
    phase = numpy.array(self._collect_fields_block(phase_file_name, phase_field_name, phase_field_components))
    self.real_part = amp*numpy.cos(phase*_deg2rad)
    self.imag_part = -amp*numpy.sin(phase*_deg2rad)

  def load_complex_components_block(self, 
      real_file_name, real_field_name, real_field_components,
      imag_file_name, imag_field_name, imag_field_components):

    self.real_part = numpy.array(self._collect_fields_block(real_file_name, real_field_name, real_field_components))
    self.imag_part = numpy.array(self._collect_fields_block(imag_file_name, imag_field_name, imag_field_components))

  def _collect_fields_block(self, file_name, field_name, field_components):
    if file_name==self.grid_file_name:
      nci = self.nci
    else:
      # copies grid,mask and ranges information from self.nci (the "grid" netCDF file)
      nci = netcdf_reader.NetCDFInterpolator(file_name, self.nci)
    nci.set_field(field_name)

    val=[]
    for component in field_components:
      val.append(nci.val[component,:,:])
    return val

  def set_time(self, t):
    if not hasattr(self,"real_part"):
      raise Exception("Need to call load_amplitudes_and_phases() first!")
    val = self.tide.from_complex_components(self.real_part, self.imag_part, t)
    self.interpolator = netcdf_reader.Interpolator(self.nci.origin, self.nci.delta, val, self.nci.mask)

  def get_val(self, x):
    if not hasattr(self,"interpolator"):
      raise Exception("Need to call set_time() first!")
    return self.interpolator.get_val(x)

def AMCGTidalInterpolator(tide, netcdf_file_name, ranges=None):
  tnci = TidalNetCDFInterpolator(tide, netcdf_file_name,
        ('latitude','longitude'),('latitude','longitude'),
        ranges=ranges)
  if "mask" in tnci.nci.nc.variables:
    tnci.set_mask("mask")

  amplitude_field_names = []
  phase_field_names = []
  for constituent in tide.constituents:
    amplitude_field_names.append( constituent.lower()+'amp' )
    phase_field_names.append( constituent.lower()+'phase' )
  tnci.load_amplitudes_and_phases(netcdf_file_name, amplitude_field_names,
      netcdf_file_name, phase_field_names)
  return tnci

def OTPSncTidalInterpolator(tide, grid_file_name, data_file_name,
    ranges=None):
  # read grid, ranges and mask from grid netCDF
  tnci = TidalNetCDFInterpolator(tide, grid_file_name,
      ('nx','ny'), ('lon_z','lat_z'), ranges=ranges)
  if "mz" in tnci.nci.nc.variables:
    tnci.set_mask("mz")
  # now swap its nci (keeping all above information) with one for the data file
  tnci.nci = netcdf_reader.NetCDFInterpolator(data_file_name, tnci.nci)

  # constituents available in the netCDF file
  constituents = tnci.nci.nc.variables['con'][:]
  # dict that maps constituent names to indices
  constituent_index = dict(((constituent.tostring().strip(' \x00').lower(),i) for i,constituent in enumerate(constituents)))
  # the indices of the requested constituents
  components = [constituent_index[constituent.lower()] for constituent in tide.constituents]
  tnci.load_complex_components_block(data_file_name, 'hRe', components,
      data_file_name, 'hIm', components)
  return tnci

  
def FESTidalInterpolator(tide, fes_file_name, ranges=None):
  # read grid, ranges and mask from grid netCDF
  tnci = TidalNetCDFInterpolator(tide, fes_file_name,
      ('Y','X'), ('lat','lon'), ranges=ranges)
  fill_value = tnci.nci.nc.variables['Ha'].missing_value
  tnci.set_mask_from_fill_value('Ha', fill_value)

  # constituents available in the netCDF file
  constituents = tnci.nci.nc.variables['spectrum'][:]
  # dict that maps constituent names to indices
  constituent_index = dict(((constituent.tostring().strip(' \x00').lower(),i) for i,constituent in enumerate(constituents)))
  # the indices of the requested constituents
  components = [constituent_index[constituent.lower()] for constituent in tide.constituents]
  tnci.load_amplitudes_and_phases_block(fes_file_name, 'Ha', components,
      fes_file_name, 'Hg', components)
  return tnci
