from __future__ import print_function
import uptide
import datetime
import numpy
from numpy.testing import assert_almost_equal, assert_equal
import math
import pytest
import pytz

pacific = pytz.timezone('US/Pacific')


# without timezone the assumption is UTC, so the following 2
# date-times should be the same
@pytest.fixture(scope="module",
                params=[datetime.datetime(2003, 1, 17, 20, 30),
                        pacific.localize(datetime.datetime(2003, 1, 17, 12, 30))])
def tide(request):
    td = uptide.Tides()
    td.set_initial_time(request.param)
    return td


@pytest.fixture(scope="module")
def mntide():
    td = uptide.Tides(['M'+str(n) for n in range(2, 13)])
    td.set_initial_time(datetime.datetime(2003, 1, 17, 19, 30))
    return td


def test_compute_nodal_corrections(tide):
    # this is a bit of a silly test: tests that no nodal correction changes
    # needs updating every time a constituent is added
    # at which point we assume that what is currently computed is correct
    tide.compute_nodal_corrections(10000)
    # test phase corrections
    assert_almost_equal(tide.u,
                        [0., -0.14205223, 0.17237799, 0.17237799, 0., 0.,
                         -0.20653437, -0.01675897, -0.03351794, 0., -0.03351794, -0.28250837,
                         -0.03351794, -0.03351794, -0.03351794, 0., -0.03351794, -0.03351794,
                         -0.03351794, 0., -0.27898729, 0., -0.06703588, -0.10055381,
                         0.03351794, -0.31602631, -0.17557017, 0.13886005, -0.03351794, -0.06703588,
                         -0.06703588, 0., -0.34954425, -0.06703588, -0.37827392, -0.03351794,
                         0., -0.37827392, -0.37827392, -0.37827392, -0.03351794, -0.37827392,
                         0., 0., -0.05027691, -0.06703588, -0.08379485, -0.10055381,
                         -0.11731278, -0.13407175, -0.15083072, -0.16758969, -0.18434866, -0.20110763])
    # test amplitude corrections
    print(tide.f)
    assert_almost_equal(tide.f,
                        [1., 1.05252903, 1.08466026, 1.08466026, 1., 1.,
                         1.06437745, 0.99251489, 0.98502979, 1., 0.98502979, 1.13971568,
                         0.98502979, 0.98502979, 0.98502979, 1., 0.98502979, 0.98502979,
                         0.98502979, 1., 1.25475539, 0.97005958, 0.97005958, 0.95508937,
                         0.98502979, 1.12474547, 1.05555882, 1.06955531, 0.98502979, 0.97005958,
                         0.97005958, 1., 1.10905669, 0.97005958, 1.21050452, 0.98502979,
                         0.94740196, 1.21050452, 1.21050452, 1.21050452, 0.93243175, 1.21050452,
                         1., 1., 0.97754468, 0.97005958, 0.96257447, 0.95508937,
                         0.94760426, 0.94011916, 0.93263405, 0.92514895, 0.91766384, 0.91017873])


def test_compute_mn_nodal_corrections(mntide):
    # test Mn nodal corrections:
    mntide.compute_nodal_corrections(1234.)
    for i in range(1, len(mntide.f)):
        n = i+2
        assert_almost_equal(mntide.f[i]-1.0, (mntide.f[0]-1.0)*(n/2.))
        assert_almost_equal(mntide.u[i], mntide.u[0]*(n/2.))


def test_ap_vs_complex(tide):
    N = len(tide.constituents)
    a = numpy.random.random_sample((N, 100))
    p = numpy.random.random_sample((N, 100))*2*math.pi
    t = 12345.
    from_ap = tide.from_amplitude_phase(a, p, t)
    real = a*numpy.cos(p)
    imag = -a*numpy.sin(p)
    from_complex = tide.from_complex_components(real, imag, t)
    assert_almost_equal(numpy.linalg.norm(from_ap-from_complex), 0.0)


def test_all_constituents_are_tested(tide, mntide):
    supported = set(uptide.tidal.omega.keys())
    tested = set(tide.constituents)
    tested.update(set(mntide.constituents))
    assert_equal(tested, supported)


def test_minimum_Rayleigh_period():
    tide = uptide.Tides(['M2', 'S2'])
    ind1, ind2 = tide.get_closest_constituents()
    assert tide.constituents[ind1] == 'M2' and tide.constituents[ind2] == 'S2'
    assert_equal(tide.get_minimum_Rayleigh_period()/86400., 14.76536317823275)

    tide = uptide.Tides(['M2', 'S2', 'N2', 'O1', 'P1', 'Q1', 'M4'])
    ind1, ind2 = tide.get_closest_constituents()
    assert tide.constituents[ind1] == 'Q1' and tide.constituents[ind2] == 'O1'
    assert_equal(tide.get_minimum_Rayleigh_period()/86400., 27.554605699504986)

    tide = uptide.Tides(['M2', 'S2', 'N2', 'K2', 'O1', 'P1', 'Q1', 'K1', 'M4', 'S1', 'MU2', 'NU2', 'L2', 'T2', 'Z0'])
    ind1, ind2 = tide.get_closest_constituents()
    assert tide.constituents[ind1] == 'T2' and tide.constituents[ind2] == 'S2'
    assert_equal(tide.get_minimum_Rayleigh_period()/86400., 365.2596411154916)


def test_select_constituents():
    constituents = ['M2', 'S2', 'N2', 'K2', 'O1', 'P1', 'Q1', 'K1', 'M4', 'S1', 'MU2', 'NU2', 'L2', 'T2', 'Z0']
    assert uptide.select_constituents(constituents, 15*86400.) == ['M2', 'S2', 'O1', 'P1', 'M4', 'Z0']
    assert uptide.select_constituents(constituents, 31*86400.) == ['M2', 'S2', 'N2', 'O1', 'P1', 'Q1', 'M4', 'Z0']
