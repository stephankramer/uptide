import unittest
import uptide.ellipse as ue
import math

class TestEllipse(unittest.TestCase):
  def test_tidal_ellipse_parameters(self):
    # first let's stick in some random numbers:
    a,b,theta = ue.tidal_ellipse_parameters(1.0, 2.0, 3.0, 4.0)
    for x,y in zip((a,b,theta),(3.0315505539181253, 0.8998340063804587, -1.4195315055223143)):
      self.assertAlmostEqual(x,y)
    # if u and v are in phase, the ellipse should be flat
    # and in the direction of the combined components
    a,b,theta = ue.tidal_ellipse_parameters(3.0,math.pi/2,4.0,math.pi/2)
    for x,y in zip((a,b,theta),(5.0, 0.0, math.atan2(4.0,3.0))):
      self.assertAlmostEqual(x,y)
    # if u and v are in anti-phase, the ellipse should have au and av
    # as its minor and major radii and the ellipse should be cartesian aligned
    a,b,theta = ue.tidal_ellipse_parameters(1.0, 1.2, 2.0, 1.2+math.pi/2.)
    for x,y in zip((a,b,math.cos(theta)),(2.0,1.0,0.0)):
      self.assertAlmostEqual(x,y)
    a,b,theta = ue.tidal_ellipse_parameters(7.0, -math.pi*2/3., 1.0, -math.pi/6.)
    for x,y in zip((a,b,math.sin(theta)),(7.0,1.0,0.0)):
      self.assertAlmostEqual(x,y)
