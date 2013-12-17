import unittest
import numpy
import uptide.constituent_combinations as uc
import uptide.ukho_nodal_corrections

class TestConstituentCombinations(unittest.TestCase):
  """Tests the uptide.constituent_combinations module"""
  def setUp(self):
    data = numpy.loadtxt('tests/ukho.csv', delimiter=',',dtype=str)
    self.omegas = dict()
    for c,omg in data[:,0:2]:
      self.omegas[c.upper()] = float(omg)

  def test_decompose_constituents(self):
    """Run decompose_constituents on all combined constituents of the ukho database,
    and check if the frequencies add up. Then run each of the constituents through
    ukho_nodal_corrections.nodal_corrections() to ensure that the decomposition can be used
    to compute nodal corrections."""

    for key,value in self.omegas.iteritems():
      if key in uc.ukho_unknown_combinations:
        self.assertRaises(uc.ConstituentCombinationException, uc.decompose_constituents, key.upper())
      else:
        letters, periods, multiplicities, diurnal = uc.decompose_constituents(key.upper())
        if len(periods)==1:
          pass
        else:
          omega = 0.0
          if not 0 in periods and len(periods)>1:
            for l,p,m in zip(letters,periods,multiplicities):
              omega += self.omegas[l+str(p)]*m
          self.assertAlmostEqual(value, omega, places=4)

          for l,p in zip(letters,periods):
             uptide.ukho_nodal_corrections.nodal_corrections_single_constituent(l+str(p), 0., 0., 0.)

if __name__ == "__main__":
    unittest.main()
