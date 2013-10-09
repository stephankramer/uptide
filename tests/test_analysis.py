import unittest
import uptide
import uptide.analysis as ua
import datetime
import numpy
import math
VISUAL = False
if VISUAL:
  from pylab import plot,show

class TestAnalysis(unittest.TestCase):
  """Tests the uptide.analysis class"""
  def setUp(self):
    self.tide = uptide.Tides()
    self.tide.set_initial_time(datetime.datetime(2003,1,17,19,30))

  def test_harmonic_analysis(self):
    N = len(self.tide.constituents)
    a = numpy.random.random_sample(N)
    p = numpy.random.random_sample(N)*2*math.pi
    trange = numpy.arange(0,86400*30,600)
    x = numpy.array([self.tide.from_amplitude_phase(a, p, t) for t in trange])
    a2, p2 = ua.harmonic_analysis(self.tide, x, trange)
    y = numpy.array([self.tide.from_amplitude_phase(a2, p2, t) for t in trange])
    if VISUAL:
      plot(trange, x)
      plot(trange, y)
      show()
    self.assertAlmostEqual(numpy.linalg.norm(x-y), 0.0)

if __name__ == '__main__':
      unittest.main()
