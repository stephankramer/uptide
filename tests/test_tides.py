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
        ['Q1', 'P1', 'K2', 'Mm', 'S2',
          'MS4', 'MN4', 'K1', 'Mf',
          'L2', 'M2', 'N2', 'Z0', 'Sa',
          'O1', 'S1', 'Ssa', '2N2', 'MU2',
          'NU2', 'T2', 'J1',  'EPS2', 'LAMBDA2',
          'MSQM', 'MFM', 'MTM', 'MSF',
          'MKS2', 'R2', 'S4', 'N4'
          ])
    self.tide.set_initial_time(datetime.datetime(2003,1,17,19,30))
    self.mntide = uptide.Tides(['M'+str(n) for n in range(2,13)])
    self.mntide.set_initial_time(datetime.datetime(2003,1,17,19,30))

  def test_compute_nodal_corrections(self):
    self.tide.compute_nodal_corrections(10000)
    # test phase corrections
    for x,y in zip(self.tide.u, 
        [ 0.17238093,  0.        , -0.28251319,  0.        ,  0.        ,
         -0.03351851, -0.06703702, -0.14205465, -0.37828036,
         -0.03351851, -0.03351851, -0.03351851,  0.        ,  0.        ,
          0.17238093,  0.        ,  0.0       , -0.03351851, -0.03351851,
         -0.03351851,  0.        , -0.20653789, -0.03351851, -0.03351851,
         -0.37828036, -0.37828036, -0.37828036, -0.03351851, 
         -0.31603169,  0.        ,  0.        , -0.06703702]):
      self.assertAlmostEqual(x,y)
    # test amplitude corrections
    for x,y in zip(self.tide.f,
        [ 1.08465367,  1.        ,  1.13970561,  0.94740654,  1.        ,
          0.98503109,  0.97006218,  1.05252498,  1.21048994,
          0.98503109,  0.98503109,  0.98503109,  1.        ,  1.        ,
          1.08465367,  1.        ,  1.        ,  0.98503109,  0.98503109,
          0.98503109,  1.        ,  1.06437150,  0.98503109,  0.98503109,
          1.21048994,  1.21048994,  1.21048994,  0.98503109,
          1.12473670,  1.        ,  1.        ,  0.97006218]):
      self.assertAlmostEqual(x,y)
    # test Mn nodal corrections:
    self.mntide.compute_nodal_corrections(1234.)
    for i in range(1,len(self.mntide.f)):
      n=i+2
      self.assertAlmostEqual(self.mntide.f[i]-1.0, (self.mntide.f[0]-1.0)*(n/2.))
      self.assertAlmostEqual(self.mntide.u[i], self.mntide.u[0]*(n/2.))

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

  def test_all_constituents_are_tested(self):
    supported = set(uptide.tidal.omega.keys())
    tested = set(self.tide.constituents)
    tested.update(set(self.mntide.constituents))
    self.assertEqual(supported, tested)

if __name__ == '__main__':
      unittest.main()
