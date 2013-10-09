import unittest
import uptide.tidal as ut
import datetime
import math

class TestTidal(unittest.TestCase):
  """Tests the uptide.tidal module (basic astronomical calculations)"""
  def setUp(self):
    # a random, insignificant date
    self.date = datetime.datetime(1975,10,24,1,2,3)
    # I can see great things in the stars for that date
    self.HshpNpp = ut.astronomical_argument(self.date)

  def test_astronomical_argument(self):
     for x,y in zip( self.HshpNpp, 
         (15.5125,
           365118.84020756715,
           27571.821754149107,
           3419.050174721317,
           -1207.0853423332171,
           282.52421802604846) ):
      self.assertAlmostEqual(x,y)

  def test_nodal_arguments(self):
    H,s,h,p,N,pp = self.HshpNpp
    # pick a lunar, solar and overtide
    f, u = ut.nodal_corrections(['M2', 'S2', 'MS4'], N, pp)
    for x,y in zip(f, ([1.0223111452791003, 1.0, 1.0223111452791003])):
      self.assertAlmostEqual(x,y)
    for x,y in zip(u, [0.02923862937459855, -0.0, 0.02923862937459855]):
      self.assertAlmostEqual(x,y)

  def test_tidal_arguments(self):
    H,s,h,p,N,pp = self.HshpNpp
    # pick a lunar, solar and overtide
    phi = ut.tidal_arguments(['N2', 'S2', 'M4'], self.date)
    for x,y in zip(phi, [-18094.924456049303, 0.54148840043124069, -23564.144432407946]):
      self.assertAlmostEqual(x,y)

  def test_omegas(self):
    self.assertAlmostEqual(ut.omega['M2'], 0.0001405189)
    # exactly twice a day:
    self.assertAlmostEqual(ut.omega['S2'], 4*math.pi/86400.)
    self.assertAlmostEqual(ut.omega['Ssa'], 3.9821275945895842e-07)

if __name__ == '__main__':
      unittest.main()
