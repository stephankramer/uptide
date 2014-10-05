from distutils.core import setup
setup(name='uptide',
      version='0.1',
      author='Stephan Kramer',
      author_email='s.kramer@imperial.ac.uk',
      description="""uptide is a python package for tidal calculations. It computes the tidal
free surface height or velocities from the amplitudes and phases of the tidal
constituents. These amplitudes and phases can be read from global tidal
solutions such as OSU/OTIS or FES2012. Some
limited functionality for tidal harmonic analysis is also available.""",
      url='http://www.opentidalfarm.com',
      packages = ['uptide'],
      classifiers = [
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)',
        'Programming Language :: Python :: 2.7',
        'Topic :: Scientific/Engineering']
      )
