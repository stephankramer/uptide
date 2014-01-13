import unittest
import os
try:
  import uptide.fes2012
except ImportError:
  pass

try:
  fes_data = os.environ['FES_DATA']
except:
  fes_data = None

if fes_data is not None:
  class TestFes(unittest.TestCase):
    def setUp(self):
      self.fes = uptide.fes2012.Fes(os.path.join(fes_data,'fes.ini'))

    def test_get_tide(self):
      h = self.fes.get_tide(58., -6., 0.)
      self.assertAlmostEqual(h, -96.5907077812)

    def test_fes_error(self):
      # probing a point in land should raise a FesException
      self.assertRaises(uptide.fes2012.FesException, self.fes.get_tide, 59, 7., 0.)
      # with the following error message
      self.assertEqual(self.fes.error(), 'Tide is undefined')

    def tearDown(self):
      # make sure the __dealloc__ method works
      del self.fes

if __name__ == '__main__':
      unittest.main()
  
