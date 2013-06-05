import math

def compute_focus_squared(au,pu,av,pv):
  """Computes the square of the focus of the ellipse from the
  amplitide and phase of x and y component of velocity. The focus f
  is the distance between the centre and each of the foci of the 
  ellipse. If a and b are the major and minor radii, then this returns:

     f**2 = a**2 - b**2"""
  return math.sqrt(au**4+av**4+2*au**2*av**2*math.cos(2*(pu-pv)))

def tidal_ellipse_parameters(au,pu,av,pv):
  """Computes the length, width and radius of a tidal ellipse
  from the amplitude and phase of the x and y components of velocity."""
  f2 = compute_focus_squared(au,pu,av,pv)
  a = math.sqrt(0.5*(au**2+av**2+f2))
  b = math.sqrt(0.5*(au**2+av**2-f2))
  # formula A3:4a from Pugh
  delta = 0.5*math.atan2(av**2*math.sin(2*(pu-pv)),au**2+av**2*math.cos(2*(pu-pv)))
  # formula A3:8 from Pugh
  theta = math.atan2(av*math.cos(pu-pv-delta), au*math.cos(delta))
  return a, b, theta
