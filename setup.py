from distutils.core import setup
from distutils.extension import Extension
import os
import sys

ext_modules = []
script_args = sys.argv[1:]
if '--with-fes' in script_args:
  script_args.remove('--with-fes')
  library_dirs = []; include_dirs = []
  if 'FES_DIR' in os.environ:
    library_dirs.append(os.path.join(os.environ['FES_DIR'], 'src'))
    include_dirs.append(os.path.join(os.environ['FES_DIR'], 'src'))
  if 'FES_LIB_DIR' in os.environ:
    library_dirs.append(os.environ['FES_LIB_DIR'])
  if 'FES_INCLUDE_DIR' in os.environ:
    include_dirs.append(os.environ['FES_INCLUDE_DIR'])
  ext_modules.append(Extension("fes2012", ["uptide/fes2012.c"],
    libraries=["fes", "netcdf"],
    library_dirs=library_dirs,
    include_dirs=include_dirs))

setup(name='uptide',
      version='0.1',
      author='Stephan Kramer',
      author_email='s.kramer@imperial.ac.uk',
      url='http://www.opentidalfarm.com',
      packages = ['uptide'],
      script_args = script_args,
      ext_package = 'uptide',
      ext_modules = ext_modules
      )
