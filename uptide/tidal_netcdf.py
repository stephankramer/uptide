import numpy
import uptide.netcdf_reader as netcdf_reader
import itertools
import os.path
import six

_deg2rad = numpy.pi/180.


class TidalNetCDFInterpolator(object):
    def __init__(self, tide, grid_file_name, dimensions, coordinate_fields,
                 ranges=None, mask=None):
        """Initiate a TidalNetCDFInterpolator. The specification of the names of the dimensions
        and coordinate_fields is the same as for the NetCDFInterpolator class, see its documentation.
        ranges and mask may be specified in a similar way to the NetCDFInterpolator class.
        NOTE: setting a correct coordinate ranges is strongly recommended when reading
        from a NetCDF data base that is significantly bigger than the region of interest,
        as otherwise the tidal signal will be reconstructed for all points of the NetCDF grid
        leading to very slow computations. In parallel the ranges should be set to the coordinate
        range of the local process domain.

        The tides may be stored as amplitudes and phases or as the real and complex components.
        For both cases they may be stored in separate fields for each constituent, or in one 3D
        field, where the third dimension corresponds to the different constituents. This third dimension
        should be the first (outer) dimension of the field as stored in the NetCDF file. The methods
        corresponding to these four storage formats are:

            tnci.load_amplitudes_and_phases("amplitude_file.nc", ("M2amp", "S2amp", "K1amp"),
                  "phase_file.nc", ("M2phase", "S2phase", "K1phase"))
            tnci.load_complex_components("realcomp_file.nc", ("M2real", "S2real", "K1real"),
                  "imagcomp_file.nc", ("M2imag", "S2imag", "K1imag"))
            tnci.load_amplitudes_and_phases_block("amplitude_file.nc", "amplitude_field", (0, 1, 2),
                  "phase_file.nc", "phase_field", (0, 1, 2))
            tnci.load_complex_components("realcomp_file.nc", "real_component", (0, 1, 2),
                  "imag_component", (0, 1, 2))

        The specified consituents should corresponds to tide.constituents. In the last two cases,
        "amplitude_field", "phase_field", "real_component", and "imag_component" refer to the
        3D field, and (0, 1, 2) which indices in the first dimension of that field corresponds to which of
        the constituents. In the first two cases instead of specifying a single file for all consituents, an array
        of filenames may also be provided to specify a separate file for each consituent.

        Three helper functions are available to create a TidalNetCDFInterpolator based on three
        different conventions for storage format: TPXOTidalInterpolator, FESTidalInterpolator and
        AMCGTidalInterpolator. These inititate a TidalNetCDFInterpolator and call the correct load_...() method.

        After this, the calling sequence is:

            tnci.set_time(t) # t in seconds after the datetime set with tide.set_initial_time()
            tnci.get_val(x)  # interpolate the tidal signal in location x

        Note that each call to set_time() the tidal signal is reconstructed in all points of the (restricted)
        NetCDF grid. Therefore this method is only efficient if a significant number of interpolations are
        done for each time.

        """
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
        """Load amplitude and phases of the different constituents where amplitudes
        and phases are stored as separate fields in the NetCDF. The amplitude and phase
        field names should be in the same order as tide.constituents. ampltide and
        phase file_name may be a single string, or an array of strings to indicate
        seperate filenames for each constituent."""

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
        """Load real and imaginary component of the different constituents where the complex
        components are stored as separate fields in the NetCDF. The complex component
        field names should be in the same order as tide.constituents.
        The file_name may be a single string, or an array of strings to indicate
        seperate filenames for each constituent."""
        self.real_part = numpy.array(self._collect_fields_val(real_file_name, real_field_names))
        self.imag_part = numpy.array(self._collect_fields_val(imag_file_name, imag_field_names))

    def _collect_fields_val(self, file_name, field_names):
        val = []
        if isinstance(file_name, six.string_types):
            file_names = itertools.repeat(file_name)
        else:
            file_names = file_name

        nci_filenm = self.grid_file_name
        nci = self.nci

        for filenm, fieldnm in zip(file_names, field_names):
            if not filenm == nci_filenm:
                # copies grid, mask and ranges information from self.nci (the "grid" netCDF file)
                nci = netcdf_reader.NetCDFInterpolator(filenm, self.nci)
                nci_filenm = filenm
            nci.set_field(fieldnm)
            if nci.dim_order[0] == 0:
                val.append(nci.val[:])
            else:
                val.append(nci.val[:].T)
        return val

    def load_amplitudes_and_phases_block(self,
                                         amplitude_file_name, amplitude_field_name, amplitude_field_components,
                                         phase_file_name, phase_field_name, phase_field_components):
        """Load amplitude and phases of the different constituents where amplitudes
        and phases are stored in a single 3d field. The first dimension of this
        field should correspond the different constituents stored. ..._field_components
        refers to which indices of this first dimension correspond to the constituents
        specified in tide.consituents."""

        amp = numpy.array(self._collect_fields_block(amplitude_file_name, amplitude_field_name, amplitude_field_components))
        phase = numpy.array(self._collect_fields_block(phase_file_name, phase_field_name, phase_field_components))
        self.real_part = amp*numpy.cos(phase*_deg2rad)
        self.imag_part = -amp*numpy.sin(phase*_deg2rad)

    def load_complex_components_block(self,
                                      real_file_name, real_field_name, real_field_components,
                                      imag_file_name, imag_field_name, imag_field_components):
        """Load complex components of the different constituents where the components
        are stored in a single 3d field. The first dimension of this
        field should correspond the different constituents stored. ..._field_components
        refers to which indices of this first dimension correspond to the constituents
        specified in tide.consituents."""

        self.real_part = numpy.array(self._collect_fields_block(real_file_name, real_field_name, real_field_components))
        self.imag_part = numpy.array(self._collect_fields_block(imag_file_name, imag_field_name, imag_field_components))

    def _collect_fields_block(self, file_name, field_name, field_components):
        if file_name == self.grid_file_name:
            nci = self.nci
        else:
            # copies grid, mask and ranges information from self.nci (the "grid" netCDF file)
            nci = netcdf_reader.NetCDFInterpolator(file_name, self.nci)
        nci.set_field(field_name)

        val = []
        for component in field_components:
            if nci.dim_order[0] == 0:
                val.append(nci.val[component, :, :])
            else:
                val.append(nci.val[component, :, :].T)
        return val

    def set_time(self, t):
        """Set the time in seconds after the datetime specified by tide.set_initial_time(). Recomputes
        the tidal signal on all points of the NetCDF grid."""
        if not hasattr(self, "real_part"):
            raise Exception("Need to call load_amplitudes_and_phases() first!")
        val = self.tide.from_complex_components(self.real_part, self.imag_part, t)
        self.interpolator = netcdf_reader.Interpolator(self.nci.origin, self.nci.delta, val, self.nci.mask)

    def get_val(self, x, allow_extrapolation=False):
        """Interpolates the tidal signal in point x, computed in set_time(). The order
        of the coordinates x is determined by the storage order in the NetCDF file."""
        if not hasattr(self, "interpolator"):
            raise Exception("Need to call set_time() first!")
        return self.interpolator.get_val(x, allow_extrapolation)


def AMCGTidalInterpolator(tide, netcdf_file_name, ranges=None):
    tnci = TidalNetCDFInterpolator(tide, netcdf_file_name,
                                   ('latitude', 'longitude'), ('latitude', 'longitude'),
                                   ranges=ranges)
    """Create a TidalNetCDFInterpolator based on the 'AMCG' storage conventions
    where amplitudes and phases are stored in separate fields in a single file
    with field names such as M2amp, M2phase, etc. If present a field named "mask"
    with 0.0 for land and 1.0 for sea will also be recognized."""
    if "mask" in tnci.nci.nc.variables:
        tnci.set_mask("mask")

    amplitude_field_names = []
    phase_field_names = []
    for constituent in tide.constituents:
        amplitude_field_names.append(constituent.lower()+'amp')
        phase_field_names.append(constituent.lower()+'phase')
    tnci.load_amplitudes_and_phases(netcdf_file_name, amplitude_field_names,
                                    netcdf_file_name, phase_field_names)
    return tnci


def TPXOTidalInterpolator(tide, grid_file_name, data_file_name,
                          ranges=None):
    """Create a TidalNetCDFInterpolator from OTPSnc NetCDF files, where
    the grid is stored in a separate file (with "lon_z", "lat_z" and "mz"
    fields). The actual data is read from a seperate file with hRe and hIm
    fields."""
    # read grid, ranges and mask from grid netCDF
    tnci = TidalNetCDFInterpolator(tide, grid_file_name,
                                   ('nx', 'ny'), ('lon_z', 'lat_z'), ranges=ranges)
    if "mz" in tnci.nci.nc.variables:
        tnci.set_mask("mz")
    # now swap its nci (keeping all above information) with one for the data file
    tnci.nci = netcdf_reader.NetCDFInterpolator(data_file_name, tnci.nci)

    # constituents available in the netCDF file
    constituents = tnci.nci.nc.variables['con'][:]
    # dict that maps constituent names to indices
    constituent_index = dict(((constituent.tostring().decode().strip(' \x00').lower(), i) for i, constituent in enumerate(constituents)))
    # the indices of the requested constituents
    components = [constituent_index[constituent.lower()] for constituent in tide.constituents]
    tnci.load_complex_components_block(data_file_name, 'hRe', components,
                                       data_file_name, 'hIm', components)
    return tnci


# old name:
OTPSncTidalInterpolator = TPXOTidalInterpolator


def FESTidalInterpolator(tide, fes_file_name, ranges=None):
    # read grid, ranges and mask from grid netCDF
    """Create a TidalNetCDFInterpolator from FES NetCDF files, where
    all constituents are stored in a single file. The amplitudes
    and phases are read from its Ha and Hg fields."""
    tnci = TidalNetCDFInterpolator(tide, fes_file_name,
                                   ('Y', 'X'), ('lat', 'lon'), ranges=ranges)
    fill_value = tnci.nci.nc.variables['Ha'].missing_value
    tnci.set_mask_from_fill_value('Ha', fill_value)

    # constituents available in the netCDF file
    constituents = tnci.nci.nc.variables['spectrum'][:]
    # dict that maps constituent names to indices
    constituent_index = dict(((constituent.tostring().decode("utf-8").strip(' \x00').lower(), i) for i, constituent in enumerate(constituents)))
    # the indices of the requested constituents
    components = [constituent_index[constituent.lower()] for constituent in tide.constituents]
    tnci.load_amplitudes_and_phases_block(fes_file_name, 'Ha', components,
                                          fes_file_name, 'Hg', components)
    return tnci


def FES2012TidalInterpolator(tide, fes_ini_file_name, fes_data_path=None, ranges=None):
    if fes_data_path is None:
        fes_data_path, tail = os.path.split(fes_ini_file_name)
    # remove double and trailing /s, change '' to '.':
    fes_data_path = os.path.normpath(fes_data_path)

    ini = read_fes_ini_file(fes_ini_file_name, fes_data_path)
    first_entry = next(ini['TIDE'].itervalues())
    grid_file_name = first_entry['FILE']
    lon_name = first_entry['LONGITUDE']
    lat_name = first_entry['LATITUDE']
    tnci = TidalNetCDFInterpolator(tide, grid_file_name, (lon_name, lat_name),
                                   (lon_name, lat_name), ranges=ranges)

    file_names = []
    amplitude_names = []
    phase_names = []
    for constituent in tide.constituents:
        file_names.append(ini['TIDE'][constituent]['FILE'])
        amplitude_names.append(ini['TIDE'][constituent]['AMPLITUDE'])
        phase_names.append(ini['TIDE'][constituent]['PHASE'])
    tnci.load_amplitudes_and_phases(file_names, amplitude_names, file_names, phase_names)
    return tnci


def read_fes_ini_file(file_name, fes_data_path):
    ini = {}
    for line in open(file_name, 'r'):
        line = line.lstrip()
        if line == '' or line.startswith(';'):
            continue

        key, value = line.split('=')
        type, name, field = key.strip().split('_', 2)

        value = value.strip()
        if field == 'FILE':
            value = value.replace('${FES_DATA}', fes_data_path)
            value = os.path.normpath(value)

        if type not in ini:
            ini[type] = {name: {field: value}}
        elif name not in ini[type]:
            ini[type][name] = {field: value}
        else:
            ini[type][name][field] = value

    return ini
