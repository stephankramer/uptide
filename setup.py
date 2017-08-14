from distutils.core import setup
from distutils.extension import Extension
import os
import sys

ext_modules = None
script_args = sys.argv[1:]
if '--with-fes' in script_args:
  script_args.remove('--with-fes')
  library_dirs = []; include_dirs = []
  if 'FES_LIB_DIR' in os.environ:
    library_dirs.append(os.environ['FES_LIB_DIR'])
  if 'FES_INCLUDE_DIR' in os.environ:
    include_dirs.append(os.environ['FES_INCLUDE_DIR'])
  if 'FES_DIR' in os.environ:
    library_dirs.append(os.path.join(os.environ['FES_DIR'], 'src'))
    include_dirs.append(os.path.join(os.environ['FES_DIR'], 'src'))
  ext_modules = [Extension("fes2012", ["fes/fes2012.c"],
    libraries=["fes", "netcdf"],
    library_dirs=library_dirs,
    include_dirs=include_dirs)]

setup(name='uptide',
      version='0.3',
      author='Stephan Kramer',
      author_email='s.kramer@imperial.ac.uk',
      description="uptide is a python package for tidal calculations. It computes the tidal " +
"free surface height or velocities from the amplitudes and phases of the tidal " +
"constituents. These amplitudes and phases can be read from global tidal " +
"solutions such as OSU/OTIS or FES2012. Some " +
"limited functionality for tidal harmonic analysis is also available.",
      url='https://github.com/stephankramer/uptide',
      packages = ['uptide'],
      keywords = ['tides', 'tidal', 'harmonic analysis'],
      classifiers = [
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)',
        'Programming Language :: Python :: 2.7',
        'Topic :: Scientific/Engineering'],
      script_args = script_args,
      ext_package = 'uptide',
      ext_modules = ext_modules
      )
