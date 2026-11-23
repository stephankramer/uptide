import numpy
import uptide.netcdf_reader as netcdf_reader
import itertools
import os.path

_deg2rad = numpy.pi/180.


class TidalNetCDFInterpolator(object):
    def __init__(self, tide, grid_file_name, dimensions, coordinate_fields,
                 ranges=None, mask=None, coordinate_periods=None):
        """Initiate a TidalNetCDFInterpolator. The specification of the names of the dimensions
        and coordinate_fields is the same as for the NetCDFGrid class, see its documentation.
        ranges, mask and coordinate_periods may be specified in a similar way to the NetCDFGrid class.
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
        self.dimensions = dimensions
        self.coordinate_fields = coordinate_fields
        self.ranges = ranges
        self.mask = mask
        self.mask_from_fill_value = None
        self.coordinate_periods = coordinate_periods
        self._create_grid()

    def _create_grid(self):
        # NetCDFGrid is immutable, so we recreate it every time ranges or mask are changed
        self.nci = netcdf_reader.NetCDFGrid(self.grid_file_name, self.dimensions, self.coordinate_fields,
                                            ranges=self.ranges, mask=self.mask,
                                            mask_from_fill_value=self.mask_from_fill_value,
                                            coordinate_periods=self.coordinate_periods)

    def set_ranges(self, ranges):
        self.ranges = ranges
        self._create_grid()

    def set_mask(self, field_name):
        self.mask = field_name
        self.mask_from_fill_value = None
        self._create_grid()

    def set_mask_from_fill_value(self, field_name, fill_value):
        self.mask = None
        self.mask_from_fill_value = (field_name, fill_value)
        self._create_grid()

    def load_amplitudes_and_phases(self, amplitude_file_name, amplitude_field_names,
                                   phase_file_name, phase_field_names, scale=1.0):
        """Load amplitude and phases of the different constituents where amplitudes
        and phases are stored as separate fields in the NetCDF. The amplitude and phase
        field names should be in the same order as tide.constituents. ampltide and
        phase file_name may be a single string, or an array of strings to indicate
        seperate filenames for each constituent. The amplitudes are multiplied by
        scale (e.g. to convert units)."""

        amp = scale*numpy.array(self._collect_fields_val(amplitude_file_name, amplitude_field_names))
        phase = numpy.array(self._collect_fields_val(phase_file_name, phase_field_names))
        self.real_part = amp*numpy.cos(phase*_deg2rad)
        self.imag_part = -amp*numpy.sin(phase*_deg2rad)
        # NOTE: I did try several things with packing everything in a single nx x ny x nc x 2 array
        # (and making sure things are  contiguous in memory in the right order)
        # and contracting it with a nc x 2 array to compute the tides at all nx x ny points
        # but everything I tried was slower than this simple version

    def load_complex_components(self, real_file_name, real_field_names,
                                imag_file_name, imag_field_names, scale=1.0):
        """Load real and imaginary component of the different constituents where the complex
        components are stored as separate fields in the NetCDF. The complex component
        field names should be in the same order as tide.constituents.
        The file_name may be a single string, or an array of strings to indicate
        seperate filenames for each constituent. Both components are multiplied by
        scale (e.g. to convert units)."""
        self.real_part = scale*numpy.array(self._collect_fields_val(real_file_name, real_field_names))
        self.imag_part = scale*numpy.array(self._collect_fields_val(imag_file_name, imag_field_names))

    def _collect_fields_val(self, file_name, field_names):
        val = []
        if isinstance(file_name, str):
            file_names = itertools.repeat(file_name)
        else:
            file_names = file_name

        grids = {self.grid_file_name: self.nci}
        for filenm, fieldnm in zip(file_names, field_names):
            if filenm not in grids:
                # copies grid, mask and ranges information from self.nci (the "grid" netCDF file)
                grids[filenm] = netcdf_reader.NetCDFGrid(filenm, self.nci)
            val.append(grids[filenm].get_values(fieldnm))
        return val

    def load_amplitudes_and_phases_block(self,
                                         amplitude_file_name, amplitude_field_name, amplitude_field_components,
                                         phase_file_name, phase_field_name, phase_field_components,
                                         scale=1.0):
        """Load amplitude and phases of the different constituents where amplitudes
        and phases are stored in a single 3d field. The first dimension of this
        field should correspond the different constituents stored. ..._field_components
        refers to which indices of this first dimension correspond to the constituents
        specified in tide.consituents. The amplitudes are multiplied by scale."""

        amp = scale*numpy.array(self._collect_fields_block(amplitude_file_name, amplitude_field_name, amplitude_field_components))
        phase = numpy.array(self._collect_fields_block(phase_file_name, phase_field_name, phase_field_components))
        self.real_part = amp*numpy.cos(phase*_deg2rad)
        self.imag_part = -amp*numpy.sin(phase*_deg2rad)

    def load_complex_components_block(self,
                                      real_file_name, real_field_name, real_field_components,
                                      imag_file_name, imag_field_name, imag_field_components,
                                      scale=1.0):
        """Load complex components of the different constituents where the components
        are stored in a single 3d field. The first dimension of this
        field should correspond the different constituents stored. ..._field_components
        refers to which indices of this first dimension correspond to the constituents
        specified in tide.consituents. Both components are multiplied by scale."""

        self.real_part = scale*numpy.array(self._collect_fields_block(real_file_name, real_field_name, real_field_components))
        self.imag_part = scale*numpy.array(self._collect_fields_block(imag_file_name, imag_field_name, imag_field_components))

    def _collect_fields_block(self, file_name, field_name, field_components):
        if file_name == self.grid_file_name:
            grid = self.nci
        else:
            # copies grid, mask and ranges information from self.nci (the "grid" netCDF file)
            grid = netcdf_reader.NetCDFGrid(file_name, self.nci)
        return [grid.get_values(field_name, component) for component in field_components]

    def set_time(self, t):
        """Set the time in seconds after the datetime specified by tide.set_initial_time(). Recomputes
        the tidal signal on all points of the NetCDF grid."""
        if not hasattr(self, "real_part"):
            raise Exception("Need to call load_amplitudes_and_phases() first!")
        val = self.tide.from_complex_components(self.real_part, self.imag_part, t)
        self.interpolator = netcdf_reader.Interpolator(self.nci.origin, self.nci.delta, val, self.nci.mask,
                                                       self.nci.periods)

    def get_val(self, x, allow_extrapolation=False):
        """Interpolates the tidal signal in point x, computed in set_time(). The order
        of the coordinates x is determined by the storage order in the NetCDF file."""
        if not hasattr(self, "interpolator"):
            raise Exception("Need to call set_time() first!")
        return self.interpolator.get_val(x, allow_extrapolation)


def AMCGTidalInterpolator(tide, netcdf_file_name, ranges=None):
    tnci = TidalNetCDFInterpolator(tide, netcdf_file_name,
                                   ('latitude', 'longitude'), ('latitude', 'longitude'),
                                   ranges=ranges, coordinate_periods=(None, 360.))
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


def _decode_constituent_name(con):
    """Convert a char array from a netCDF 'con' variable into a lower case string."""
    return numpy.asarray(con).tobytes().decode().strip(' \x00').lower()


def _constituent_indices(con_variable, tide):
    """Return the indices in the netCDF 'con' variable (or similar) of tide.constituents."""
    constituent_index = dict((_decode_constituent_name(constituent), i) for i, constituent in enumerate(con_variable[:]))
    return [constituent_index[constituent.lower()] for constituent in tide.constituents]


# recognized units in TPXO files: maps normalized units string to (quantity, factor to convert to SI)
_tpxo_units_table = {
    'm': ('elevation', 1.0),
    'mm': ('elevation', 1e-3),
    'm/s': ('velocity', 1.0),
    'cm/s': ('velocity', 1e-2),
    'm^2/s': ('transport', 1.0),
    'cm^2/s': ('transport', 1e-4),
}


# units assumed if no units attribute is present, following the OTPSnc conventions
_tpxo_default_units = {'hRe': 'meter', 'ua': 'cm/s', 'va': 'cm/s', 'Ua': 'meter^2/s', 'Va': 'meter^2/s',
                       'URe': 'meter^2/s', 'VRe': 'meter^2/s'}


def _tpxo_units(variable, field_name):
    """Return (quantity, scale) for a TPXO netCDF variable based on its units attribute,
    where quantity is 'elevation', 'velocity' or 'transport' and scale is the factor that converts
    the values to SI units (m, m/s or m^2/s)."""
    units = getattr(variable, 'units', None) or _tpxo_default_units.get(field_name)
    if units is None:
        raise ValueError("TPXO netCDF variable {} has no units attribute".format(field_name))
    if isinstance(units, bytes):
        units = units.decode()
    normalized = units.strip().lower().replace(' ', '')
    for long_name, short_name in (('centimeter', 'cm'), ('millimeter', 'mm'), ('meter', 'm'), ('sec', 's')):
        normalized = normalized.replace(long_name, short_name)
    if normalized not in _tpxo_units_table:
        raise ValueError("Unrecognized units '{}' in TPXO netCDF file".format(units))
    return _tpxo_units_table[normalized]


def _tpxo_set_mask(tnci, grid_field_name):
    """Set the land mask from the "mX" field, or if not present (in TPXO atlas grid files)
    from the bathymetry "hX" field which is zero on land."""
    variables = tnci.nci.nc.variables
    if 'm{}'.format(grid_field_name) in variables:
        tnci.set_mask('m{}'.format(grid_field_name))
    elif 'h{}'.format(grid_field_name) in variables:
        tnci.set_mask_from_fill_value('h{}'.format(grid_field_name), 0.)


class _TPXODataFiles(object):
    """Describes the TPXO data file(s) for tide.constituents: either a single file with all constituents,
    or, if data_file_name contains '{con}', a separate file for each constituent."""
    def __init__(self, tide, data_file_name):
        self.per_constituent = '{con}' in data_file_name
        if self.per_constituent:
            self.file_names = []
            for constituent in tide.constituents:
                file_name = data_file_name.format(con=constituent.lower())
                if not os.path.exists(file_name):
                    raise FileNotFoundError("TPXO data file {} for constituent {} not found".format(file_name, constituent))
                nc = netcdf_reader.NetCDFFile(file_name, 'r')
                if 'con' in nc.variables:
                    file_constituent = _decode_constituent_name(nc.variables['con'][:])
                    if file_constituent != constituent.lower():
                        raise ValueError("TPXO data file {} contains constituent {}, expected {}".format(
                            file_name, file_constituent, constituent.lower()))
                nc.close()
                self.file_names.append(file_name)
            self.first_file_name = self.file_names[0]
        else:
            self.file_name = self.first_file_name = data_file_name
            nc = netcdf_reader.NetCDFFile(data_file_name, 'r')
            self.components = _constituent_indices(nc.variables['con'], tide)
            nc.close()

    def units(self, field_name):
        """Return (quantity, scale) of field_name, see _tpxo_units(). For per constituent files,
        these should be the same in all files."""
        file_names = self.file_names if self.per_constituent else [self.file_name]
        result = None
        for file_name in file_names:
            nc = netcdf_reader.NetCDFFile(file_name, 'r')
            units = _tpxo_units(nc.variables[field_name], field_name)
            nc.close()
            if result is None:
                result = units
            elif units != result:
                raise ValueError("Inconsistent units of {} in TPXO data files".format(field_name))
        return result

    def variables(self):
        """Names of variables available in the (first) data file."""
        nc = netcdf_reader.NetCDFFile(self.first_file_name, 'r')
        names = set(nc.variables.keys())
        nc.close()
        return names

    def load_complex_components(self, tnci, real_field_name, imag_field_name, scale):
        if self.per_constituent:
            n = len(self.file_names)
            tnci.load_complex_components(self.file_names, [real_field_name]*n,
                                         self.file_names, [imag_field_name]*n, scale=scale)
        else:
            tnci.load_complex_components_block(self.file_name, real_field_name, self.components,
                                               self.file_name, imag_field_name, self.components, scale=scale)

    def load_amplitudes_and_phases(self, tnci, amplitude_field_name, phase_field_name, scale):
        if self.per_constituent:
            n = len(self.file_names)
            tnci.load_amplitudes_and_phases(self.file_names, [amplitude_field_name]*n,
                                            self.file_names, [phase_field_name]*n, scale=scale)
        else:
            tnci.load_amplitudes_and_phases_block(self.file_name, amplitude_field_name, self.components,
                                                  self.file_name, phase_field_name, self.components, scale=scale)


def TPXOTidalInterpolator(tide, grid_file_name, data_file_name,
                          ranges=None):
    """Create a TidalNetCDFInterpolator from OTPSnc NetCDF files, where
    the grid is stored in a separate file (with "lon_z", "lat_z" and either a "mz"
    mask field, or "hz" bathymetry field with zeros on land).

    The actual data is read from hRe and hIm fields. These may either be stored in
    a single data file containing all constituents, e.g.:

        TPXOTidalInterpolator(tide, 'grid_tpxo9.nc', 'h_tpxo9.v1.nc')

    or in separate files per constituent, as in the TPXO atlas versions. In that case
    data_file_name should be a format string containing '{con}', which is substituted
    by the lower case constituent name, e.g.:

        TPXOTidalInterpolator(tide, 'grid_tpxo10atlas_v2.nc', 'h_{con}_tpxo10_atlas_30_v2.nc')

    The units of the data (meter or millimeter) are read from the data file(s), and
    returned elevations are always in meters.

    NOTE: specifying ranges is strongly recommended, in particular for the high
    resolution atlas versions."""
    # read grid, ranges and mask from grid netCDF
    tnci = TidalNetCDFInterpolator(tide, grid_file_name,
                                   ('nx', 'ny'), ('lon_z', 'lat_z'), ranges=ranges,
                                   coordinate_periods=(360., None))
    _tpxo_set_mask(tnci, 'z')

    data_files = _TPXODataFiles(tide, data_file_name)
    quantity, scale = data_files.units('hRe')
    if quantity != 'elevation':
        raise ValueError("Units of hRe field in TPXO data file should be meter or millimeter")
    data_files.load_complex_components(tnci, 'hRe', 'hIm', scale)
    return tnci


def TPXOncTidalComponentInterpolator(tide, grid_file_name, data_file_name,
                                     grid_field_name, field_name, ranges=None):
    """Create a TidalNetCDFInterpolator for a velocity or transport component from OTPSnc
    NetCDF files, where the grid is stored in a separate file (with "lon_X" and "lat_X"
    fields, and either a "mX" mask field or "hX" bathymetry field with zeros on land),
    where X=grid_field_name is u or v.

    The field_name should be either 'u' or 'v' for the eastward or northward depth averaged
    velocity, in which case values are returned in m/s, or 'U' or 'V' for the eastward or
    northward (depth integrated) transport, in which case values are returned in m^2/s.

    NOTE: in uptide versions up to 1.2, velocities were returned in cm/s (the units used
    in the TPXO files) instead of m/s.

    Similar to TPXOTidalInterpolator(), the data may either be stored in a single file
    for all constituents, or in a separate file per constituent, as in the TPXO atlas
    versions, by specifying a format string containing '{con}' as data_file_name, e.g.:

        TPXOncTidalComponentInterpolator(tide, 'grid_tpxo10atlas_v2.nc', 'u_{con}_tpxo10_atlas_30_v2.nc',
                                         'u', 'u', ranges=ranges)

    The data file(s) should contain either amplitudes and phases (e.g. "ua" and "up"), or complex
    components (e.g. "URe" and "UIm" or "uRe" and "uIm"). Whether these represent velocities or
    transports is determined from their units attribute. If depth averaged velocities are requested
    but only transports are available (as in the TPXO atlas versions), the velocities are computed by
    dividing the transport by the depth "hX" from the grid file."""
    if field_name not in ('u', 'v', 'U', 'V'):
        raise ValueError("field_name should be one of 'u', 'v', 'U' or 'V'")
    # read grid, ranges and mask from grid netCDF
    tnci = TidalNetCDFInterpolator(tide, grid_file_name,
                                   ('nx', 'ny'),
                                   ('lon_{}'.format(grid_field_name),
                                    'lat_{}'.format(grid_field_name)), ranges=ranges,
                                   coordinate_periods=(360., None))
    _tpxo_set_mask(tnci, grid_field_name)

    data_files = _TPXODataFiles(tide, data_file_name)
    variables = data_files.variables()

    # candidate (load method, field names) for both upper and lower case names - each
    # may store either velocity or transport (determined from the units)
    candidates = []
    for name in (field_name, field_name.swapcase()):
        candidates.append(('complex', name + 'Re', name + 'Im'))
        candidates.append(('ampphase', name + 'a', name.lower() + 'p'))
    available = {}
    for method, field1, field2 in candidates:
        if field1 in variables and field2 in variables:
            quantity, scale = data_files.units(field1)
            available.setdefault(quantity, (method, field1, field2, scale))

    requested = 'velocity' if field_name.islower() else 'transport'
    if requested in available:
        method, field1, field2, scale = available[requested]
        divide_by_depth = False
    elif requested == 'velocity' and 'transport' in available:
        method, field1, field2, scale = available['transport']
        divide_by_depth = True
    else:
        raise ValueError("No {} data for component {} found in TPXO data file {}".format(
            requested, field_name, data_files.first_file_name))

    if method == 'complex':
        data_files.load_complex_components(tnci, field1, field2, scale)
    else:
        data_files.load_amplitudes_and_phases(tnci, field1, field2, scale)

    if divide_by_depth:
        depth_field_name = 'h{}'.format(grid_field_name)
        if depth_field_name not in tnci.nci.nc.variables:
            raise ValueError("Need depth field {} in TPXO grid file to compute velocities from transports".format(
                depth_field_name))
        depth_units = getattr(tnci.nci.nc.variables[depth_field_name], 'units', 'meter')
        if isinstance(depth_units, bytes):
            depth_units = depth_units.decode()
        if depth_units.strip().lower() not in ('m', 'meter', 'meters'):
            raise ValueError("Unrecognized units '{}' of depth field {}".format(depth_units, depth_field_name))
        depth = numpy.ma.filled(tnci.nci.get_values(depth_field_name), 0.).astype(float)
        # land points (zero depth) are masked out anyway, avoid division by zero there
        wet = depth > 0.
        inv_depth = numpy.zeros_like(depth)
        inv_depth[wet] = 1./depth[wet]
        tnci.real_part = tnci.real_part*inv_depth
        tnci.imag_part = tnci.imag_part*inv_depth

    return tnci


# old name:
OTPSncTidalInterpolator = TPXOTidalInterpolator
OTPSncTidalComponentInterpolator = TPXOncTidalComponentInterpolator


def FESTidalInterpolator(tide, fes_file_name, ranges=None):
    # read grid, ranges and mask from grid netCDF
    """Create a TidalNetCDFInterpolator from FES NetCDF files, where
    all constituents are stored in a single file. The amplitudes
    and phases are read from its Ha and Hg fields."""
    tnci = TidalNetCDFInterpolator(tide, fes_file_name,
                                   ('Y', 'X'), ('lat', 'lon'), ranges=ranges,
                                   coordinate_periods=(None, 360.))
    fill_value = tnci.nci.nc.variables['Ha'].missing_value
    tnci.set_mask_from_fill_value('Ha', fill_value)

    # the indices of the requested constituents in the netCDF file
    components = _constituent_indices(tnci.nci.nc.variables['spectrum'], tide)
    tnci.load_amplitudes_and_phases_block(fes_file_name, 'Ha', components,
                                          fes_file_name, 'Hg', components)
    return tnci


def FES2012TidalInterpolator(tide, fes_ini_file_name, fes_data_path=None, ranges=None):
    if fes_data_path is None:
        fes_data_path, tail = os.path.split(fes_ini_file_name)
    # remove double and trailing /s, change '' to '.':
    fes_data_path = os.path.normpath(fes_data_path)

    ini = read_fes_ini_file(fes_ini_file_name, fes_data_path)
    first_entry = next(iter(ini["TIDE"].values()))
    grid_file_name = first_entry['FILE']
    lon_name = first_entry['LONGITUDE']
    lat_name = first_entry['LATITUDE']
    tnci = TidalNetCDFInterpolator(tide, grid_file_name, (lon_name, lat_name),
                                   (lon_name, lat_name), ranges=ranges,
                                   coordinate_periods=(360., None))

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
