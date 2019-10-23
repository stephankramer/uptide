"""Implements the Tides class"""
from __future__ import print_function
import uptide.tidal as tidal
import datetime
import numpy


class Tides(object):

    """Class for tidal computations."""

    def __init__(self, constituents=None):
        """Initialise tidal computation class. If no constituents are specified
        all known constituents are included. The used constituents can be queried
        from the constituents attribute."""
        if constituents is not None:
            self.constituents = [c.upper() for c in constituents]
        else:
            self.constituents = list(tidal.omega.keys())

        try:
            self.omega = numpy.array(
                [tidal.omega[constituent] for constituent in self.constituents])
        except KeyError:
            print("*** Unsupported constituent(s) in uptide: ", set(self.constituents)-set(tidal.omega.keys()))
            raise
        self.phi = None  # Greenwich argument
        self.f = None  # nodal amplitude correction
        self.u = None  # nodal phase correction

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
        H, s, h, p, N, pp = tidal.astronomical_argument(time)
        self.f, self.u = tidal.nodal_corrections(self.constituents, N, pp)

    def from_amplitude_phase(self, amplitudes, phases, t):
        """Compute the tide from provided amplitudes and phases (same order
        as constituents provided at initialisation of the object) at time t.

        The amplitudes and phases may be a list or array with
        a single value per constituent, in which case a single tidal
        value is calculated. They may also be a list of numpy arrays of any shape,
        or a numpy array whose first dimension corresponds to the number
        of consituents. In this case multiple tidal values are calculated at once."""
        # we use eta here, but this may of course just as well be velocities
        eta = 0.0
        for f, amplitude, omega, phase, phi, u in zip(self.f, amplitudes, self.omega,
                                                      phases, self.phi, self.u):
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
        for f, omega, phi, u, real_part, imag_part in zip(self.f, self.omega,
                                                          self.phi, self.u, real_parts, imag_parts):
            eta += f*(numpy.cos(omega*t+phi+u)*real_part
                      - numpy.sin(omega*t+phi+u)*imag_part)
        return eta

    def get_closest_constituents(self):
        """Return the indices of the two constituents with the closest frequency.

        Example:
        ind1, ind2 = tide.get_closest_constituents()
        print("The two closest constituents are: ", tide.constituents[ind1], tide.constituents[ind2])
        print("with periods: ", 2*pi/tide.omega[ind1], 2*pi/tide.omega[ind2])
        """
        ind = numpy.argsort(self.omega)
        diff_omega = self.omega[ind][1:] - self.omega[ind[:-1]]
        mind = diff_omega.argmin()
        return ind[mind], ind[mind+1]

    def get_minimum_Rayleigh_period(self):
        """Compute the minimum period (seconds) of time needed to distinguish the specified constituents
        according to the Rayleigh criterion."""
        ind1, ind2 = self.get_closest_constituents()
        return 2*numpy.pi/(self.omega[ind2]-self.omega[ind1])


def select_constituents(constituents, period):
    """Select constituents according to Rayleigh criterion.

    Selects those constituents such that the minimum Rayleigh period (of the two constituents
    with the closest frequency) is smaller than the specified period. Constituents should be supplied
    in order of importance. Whenever two constituents are discovered whose
    frequencies are two close the one that is further in the list is removed."""

    tide = Tides(constituents)
    min_period = tide.get_minimum_Rayleigh_period()
    if min_period < period:
        return constituents
    else:
        ind1, ind2 = tide.get_closest_constituents()
        max_ind = max(ind1, ind2)
        return select_constituents(constituents[:max_ind]+constituents[max_ind+1:], period)
