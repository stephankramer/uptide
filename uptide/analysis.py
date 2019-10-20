import numpy
import numpy.linalg


def harmonic_analysis(tide, x, t):
    """Perform tidal harmonic analysis for a given signal x at times t.
    Returns the amplitudes and phases of the constituents defined in tide
    (a Tides object), in the order of tide.constituents. The times t are
    in seconds after the date time set with tide.set_initial_time()."""
    if not len(x) == len(t):
        raise Exception("Length of x and t should be the same")
    N = len(t)
    M = len(tide.omega)

    """We first target for a solution of the form:
        eta = Z0 + \\sum_n B_n cos omega_n t + C_n sin omega_n t
                = Z0 + \\sum_n Re[ (B_n-i C_n) e^(i omega_n t) ]"""
    if "Z0" in tide.constituents:
        """we solve for the least squares approximation
                y=[B_1, B_2,...B_M, C_1, C_2,...C_M]
            that closest approximates x at all t, where for some j:
            B_j=Z0 and we leave out C_j (as its associated sin() is always zero)"""
        # the indices of the constituents other than Z0
        nonz0 = numpy.array(tide.constituents) != 'Z0'
        mat = numpy.hstack([
            numpy.cos(numpy.outer(t, tide.omega)),
            numpy.sin(numpy.outer(t, tide.omega[nonz0]))
        ])
    else:
        """we solve for the least squares approximation
                y=[Z0, B_1, B_2,...B_M, C_1, C_2,...C_M]
            that closest approximates x at all t"""
        mat = numpy.hstack([
            numpy.ones((N, 1)),
            numpy.cos(numpy.outer(t, tide.omega)),
            numpy.sin(numpy.outer(t, tide.omega))
        ])

    y = numpy.linalg.lstsq(mat, x, rcond=None)

    if "Z0" in tide.constituents:
        B = y[0][0:M]
        # the C_j associated with Z0 is left zero
        C = numpy.zeros(M)
        C[nonz0] = y[0][M:]
    else:
        B = y[0][1:M+1]
        C = y[0][M+1:]
    A = B - 1j*C
    """Now with A_n=B_n -i C_n, we have
          eta = Z0 + \\sum_n Re[ A_n e^(i omega_n t) ]
                  = Z0 + \\sum_n Re[ C_n f_n e^(i (-g+phi+u)) e^(i omega_n t)]
                  = Z0 + \\sum_n a_n f_n cos(omega_n t-g+phi_n+u_n),
    i.o.w. A_n = a_n f_n e^(i (-g+phi+u))"""
    a = numpy.abs(A)/tide.f
    arg = numpy.angle(A)
    g = (tide.phi+tide.u-arg) % (2*numpy.pi)
    return a, g


def error_analysis(mod_amp, mod_phase, obs_amp, obs_phase):
    """Perform error analysis of model and observations (or two models) based
    on amplitudes and phases of components. The test is based on
    Cummins, Patrick F., and Lie-Yauw Oey. 1997. Simulation of Barotropic
    and Baroclinic Tides off Northern British Columbia.
    Journal of Physical Oceanography 27(5). 762-81.
    The formula is roughly:
    D_n = sqrt(Ave(eta_o - eta_m)^2))
            = sqrt(0.5(h_o^2 + h_m^2) - h_o h_m cos(phi_o - phi_m))
    where h is amplitude, phi is phase. _o is observation, _m is model.
    eta is FS height, and the average is over a tidal period. D_n give
    the RMS error of the FS height.

    Input:
            tide - tide object
            mod_amp - numpy array. Model amplitudes in the order of tide.constituents
            mod_phase - numpy array. Model phases in order of tides.constiuents
            obs_amp - numpy array. Observed amplitudes in the order of tide.constituents
            obs_phase - numpy array. Observed phases in order of tide.constiuents

    Returns list of floats for D_n for each tidal constituent in tide.constiuents
    """
    # seperating this out to make it easier to read and understand...See formula above.
    D = 0.5 * ((mod_amp * mod_amp) + (obs_amp * obs_amp))
    D = D - (obs_amp*mod_amp * numpy.cos(obs_phase - mod_phase))
    D_n = numpy.sqrt(D)

    return D_n
