from setuptools import setup
import os
import sys

script_args = sys.argv[1:]

setup(name='uptide',
      version='1.0',
      author='Stephan Kramer',
      author_email='s.kramer@imperial.ac.uk',
      description="uptide is a python package for tidal calculations. It computes the tidal " +
"free surface height or velocities from the amplitudes and phases of the tidal " +
"constituents. These amplitudes and phases can be read from global tidal " +
"solutions such as TPXO  or FES2014. Some " +
"limited functionality for tidal harmonic analysis is also available.",
      url='https://github.com/stephankramer/uptide',
      packages = ['uptide'],
      keywords = ['tides', 'tidal', 'harmonic analysis'],
      classifiers = [
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering'],
      script_args = script_args,
      ext_package = 'uptide',
      )
