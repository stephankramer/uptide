"""Implements the Tides class"""
import tidal
import datetime
import numpy
import math

class Tides(object):
  """Class for tidal computations."""
  def __init__(self, constituents=None):
    """Initialise tidal computation class. If no components are specified
    all known components are included."""
    if constituents is not None:
      self.constituents = constituents
    else:
      self.constituents = tidal.omega.keys()

    self.omega = numpy.array(
        [tidal.omega[constituent] for constituent in self.constituents])
    self.phi = None # Greenwich argument
    self.f = None # nodal amplitude correction
    self.u = None # nodal phase correction

  def set_initial_time(self, datetime0):
    """Set the initial date and time by supplying its datetime object.
    In all subsequent computations the time t is specified as seconds
    since this date+time. 
    
    The (slowly varying) nodal corrections by default are 
    only computed for this date+time. If more frequent recomputation
    is required (for longer period computations), call 
    compute_nodal_corrections() at the appropriate times.

    This date+time also determines the point at which the very long term,
    (think geological timescales) non-linear (quadratic and cubic) changes
    in the normal linear angular speeds of the sun and moon orbits are
    computed. This date+time should therefore be chosen reasonably
    close (in the same era) to the times that are actually computed.""" 
    self.datetime0 = datetime0
    self.compute_nodal_corrections(0.0)
    self.phi = tidal.tidal_arguments(self.constituents, datetime0)

  def compute_nodal_corrections(self, t):
    """(Re)Compute the nodal corrections (amplitude corrections 'f'
    and phase corrections 'u') at the specified time t (in seconds
    since the date+time set with set_initial_time())"""
    time = self.datetime0 + datetime.timedelta(seconds=t)
    H,s,h,p,N,pp = tidal.astronomical_argument(time)
    self.f, self.u = tidal.nodal_corrections(self.constituents, N, pp)

  def from_amplitude_phase(self, amplitudes,phases,t):
    """Compute the tide from provided amplitudes and phases (same order
    as constituents provided at initialisation of the object) at time t.

    The amplitudes and phases may be a list or array with 
    a single value per constituent, in which case a single tidal
    value is calculated. They may also be a list of numpy arrays of any shape,
    or a numpy array whose first dimension corresponds to the number 
    of consituents. In this case multiple tidal values are calculated at once."""
    # we use eta here, but this may of course just as well be velocities
    eta = 0.0
    for f,amplitude,omega,phase,phi,u in zip(self.f, amplitudes, self.omega, 
                                       phases,self.phi, self.u):
      eta += f*amplitude*numpy.cos(omega*t-phase+phi+u)
    return eta

  def from_complex_components(self, real_parts, imag_parts, t):
    """Compute the tide from provided real and imaginary parts of the 
    tidal constituents (same order as provided at initialisation of the object) at time t.

    The real and imaginary parts may be a list or array with 
    a single value per constituent, in which case a single tidal
    value is calculated. They may also be a list of numpy arrays of any shape,
    or a numpy array whose first dimension corresponds to the number 
    of consituents. In this case multiple tidal values are calculated at once."""
    eta = 0.0
    for f,omega,phi,u,real_part,imag_part in zip(self.f, self.omega, 
                              self.phi, self.u, real_parts, imag_parts):
      eta += f*(numpy.cos(omega*t+phi+u)*real_part
        -numpy.sin(omega*t+phi+u)*imag_part)
    return eta
