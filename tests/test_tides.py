import unittest
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
        params=[datetime.datetime(2003,1,17,20,30),
         pacific.localize(datetime.datetime(2003,1,17,12,30))])
def tide(request):
    # we ask for all constituents explicitly, so that if we add new
    # constituents all tests below, except test_all_constituents_tested,
    # should still pass
    td = uptide.Tides(
        ['Q1', 'P1', 'K2', 'Mm', 'S2',
          'MS4', 'MN4', 'K1', 'Mf',
          'L2', 'M2', 'N2', 'Z0', 'Sa',
          'O1', 'S1', 'Ssa', '2N2', 'MU2',
          'NU2', 'T2', 'J1',  'EPS2', 'LAMBDA2',
          'MSQM', 'MFM', 'MTM', 'MSF',
          'MKS2', 'R2', 'S4', 'N4'
          ])
    td.set_initial_time(request.param)
    return td

@pytest.fixture(scope="module")
def mntide():
    td = uptide.Tides(['M'+str(n) for n in range(2,13)])
    td.set_initial_time(datetime.datetime(2003,1,17,19,30))
    return td

def test_compute_nodal_corrections(tide):
    tide.compute_nodal_corrections(10000)
    # test phase corrections
    for x,y in zip(tide.u, 
            [ 0.17237799,  0.        , -0.28250837,  0.        ,  0.        , -0.03351794,
             -0.06703588, -0.14205223, -0.37827392, -0.03351794, -0.03351794, -0.03351794,
              0.        ,  0.        ,  0.17237799,  0.        ,  0.        , -0.03351794,
             -0.03351794, -0.03351794,  0.        , -0.20653437, -0.03351794, -0.03351794,
             -0.37827392, -0.37827392, -0.37827392, -0.03351794, -0.31602631,  0.        ,  
              0.        , -0.06703588]):
      assert_almost_equal(x,y)
    # test amplitude corrections
    print tide.f
    for x,y in zip(tide.f,
            [ 1.08466026,  1.        ,  1.13971568,  0.94740196,  1.        ,  0.98502979, 
              0.97005958,  1.05252903,  1.21050452,  0.98502979,  0.98502979,  0.98502979,
              1.        ,  1.        ,  1.08466026,  1.        ,  1.        ,  0.98502979,
              0.98502979,  0.98502979,  1.        ,  1.06437745,  0.98502979,  0.98502979,
              1.21050452,  1.21050452,  1.21050452,  0.98502979,  1.12474547,  1.        ,
              1.        ,  0.97005958]):
      assert_almost_equal(x,y)

def test_compute_mn_nodal_corrections(mntide):
    # test Mn nodal corrections:
    mntide.compute_nodal_corrections(1234.)
    for i in range(1,len(mntide.f)):
      n=i+2
      assert_almost_equal(mntide.f[i]-1.0, (mntide.f[0]-1.0)*(n/2.))
      assert_almost_equal(mntide.u[i], mntide.u[0]*(n/2.))

def test_ap_vs_complex(tide):
    N = len(tide.constituents)
    a = numpy.random.random_sample((N,100))
    p = numpy.random.random_sample((N,100))*2*math.pi
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
