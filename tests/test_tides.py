import unittest
import uptide
import datetime
import numpy
import math

class TestTidal(unittest.TestCase):
  """Tests the uptide.Tides class"""
  def setUp(self):
    # we ask for all constituents explicitly, so that if we add new
    # constituents all tests below, except test_all_constituents_tested,
    # should still pass
    self.tide = uptide.Tides(
        ['Q1',
          'P1',
          'K2',
          'Mm',
          'S2',
          'MS4',
          'MN4',
          'M4',
          'K1',
          'Mf',
          'L2',
          'M2',
          'N2',
          'Z0',
          'Sa',
          'O1',
          'Ssa'])
    self.tide.set_initial_time(datetime.datetime(2003,1,17,19,30))

  def test_compute_nodal_corrections(self):
    self.tide.compute_nodal_corrections(10000)
    # test phase corrections
    for x,y in zip(self.tide.u,
        [ 0.17238093,  0.        , -0.28251319,  0.        ,  0.        ,
          -0.03351851, -0.06703702, -0.06703702, -0.14205465, -0.37828036,
          0.        , -0.03351851, -0.03351851,  0.        ,  0.        ,
          0.17238093,  0.        ]):
      self.assertAlmostEqual(x,y)
    # test amplitude corrections
    for x,y in zip(self.tide.f,
        [ 1.08465367,  1.        ,  1.13970561,  0.94740654,  1.        ,
          0.98503109,  0.97028625,  0.97028625,  1.05252498,  1.21048994,
          1.        ,  0.98503109,  0.98503109,  1.        ,  1.        ,
          1.08465367,  1.        ]):
      self.assertAlmostEqual(x,y)

  def test_ap_vs_complex(self):
    N = len(self.tide.constituents)
    a = numpy.random.random_sample((N,100))
    p = numpy.random.random_sample((N,100))*2*math.pi
    t = 12345.
    from_ap = self.tide.from_amplitude_phase(a, p, t)
    real = a*numpy.cos(p)
    imag = -a*numpy.sin(p)
    from_complex = self.tide.from_complex_components(real, imag, t)
    self.assertAlmostEqual(numpy.linalg.norm(from_ap-from_complex), 0.0)

