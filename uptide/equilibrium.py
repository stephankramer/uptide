import numpy
import datetime
from . import tidal


# amplitudes for equilibrium tide from table 6.1.1 in K&C
# who cite Desai '96 for these values
# NOTE: we haven't included diurnal constituents with
# amplitude smaller than J1, and not inculded R2 either (also tiny)
equilibirum_amplitudes = {
    'MF': 0.042017,
    'MM': 0.022191,
    'SSA': 0.019542,
    'MT': 0.008049,  # name in K&C, same as following 2:
    'MFM': 0.008049,  # name in FES
    'MTM': 0.008049,  # name in UKHO
    'MSM': 0.004239,
    'MSF': 0.003678,
    'SA': 0.003104,
    'K1': 0.142408,
    'O1': 0.101266,
    'P1': 0.047129,
    'Q1': 0.019387,
    'M1': 0.007965,
    'J1': 0.007965,
    'M2': 0.244102,
    'S2': 0.113572,
    'N2': 0.046735,
    'K2': 0.030875,
    'NU2': 0.008877,
    'MU2': 0.007463,
    'L2': 0.006899,
    'T2': 0.006636,
    '2N2': 0.006184,
    'EPS2': 0.001804,
    'LAMBDA2': 0.001800,
    'ETA2': 0.001727,
    # 'M3': 0.003198,  # leaving this out as I don't know the scaling of the P3 Legendre polynomial that was used
}

ALL_EQUILIBRIUM_TIDAL_CONSTITUENTS = equilibirum_amplitudes.keys()


def equilibrium_tide(tide, lat, lon, t):
    eta = 0.
    time = tide.datetime0 + datetime.timedelta(seconds=t)
    chi_list = tidal.tidal_arguments(tide.constituents, time)
    cosl2 = numpy.cos(lat) ** 2
    P = [0]*3  # Legendre polynomials
    P[0] = 3 * cosl2 / 2 - 1  # NOTE: this is half of eqn (4) in Schwiderski (Rev. of G&SP, '80)
    # but is consistent with K&C and also the original Cartwright and Tayler
    # and also what is given in:
    # Schwiderski, Global Ocean Tides. Part X. The Fortnightly Lunar Tide (Mf) Atlas of Tidal Charts and Maps.
    # 1982, page 3 - part of a series of technical US Defense reports
    # So I think Schwiderski '80 is wrong (?!)
    P[1] = numpy.sin(2*lat)
    P[2] = numpy.cos(lat)**2
    for c, chi, f, u in zip(tide.constituents, chi_list, tide.f, tide.u):
        m = int(tidal.lunar_doodson_numbers[c][0])
        eta += equilibirum_amplitudes[c] * P[m] * f * numpy.cos(tidal.omega[c]*t + m*lon + chi + u)
    return eta
