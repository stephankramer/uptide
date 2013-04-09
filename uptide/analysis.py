import numpy
import numpy.linalg

def harmonic_analysis(tide, x, t):
  if len(x)!=len(t):
    raise Exception("Length of x and t should be the same")
  N = len(t)
  M = len(tide.omega)

  # if Z0 is not indicated as one of the consituents, then we subtract
  # the average beforehand - note that this may give a slightly different
  # answer for the other constituents as well
  if not "Z0" in x:
    x = numpy.array(x)
    x = x - x.mean()

  """We first target for a solution of the form:
    eta = \sum_n B_n cos omega_n t + C_n sin omega_n t
        = \sum_n Re[ (B_n-i C_n) e^(i omega_n t) ],
  we solve for the least squares approximation
     y=[B_1,B_2,...B_M,C_1,C_2,...C_M]
  that closest approximates x at all t"""
  mat = numpy.hstack([
    numpy.cos(numpy.outer(t,tide.omega)), 
    numpy.sin(numpy.outer(t,tide.omega))
    ])

  y = numpy.linalg.lstsq(mat, x)

  B = y[0][0:M]
  C = y[0][M:]
  A = B - 1j*C
  """Now with A_n=B_n -i C_n, we have
     eta = \sum_n Re[ A_n e^(i omega_n t) ]
         = \sum_n Re[ C_n f_n e^(i (-g+phi+u)) e^(i omega_n t)]
         = \sum_n a_n f_n cos(omega_n t-g+phi_n+u_n),
  i.o.w. A_n = a_n f_n e^(i (-g+phi+u))"""
  a = numpy.abs(A)/tide.f
  arg = numpy.angle(A)
  g = tide.phi+tide.u-arg
  return a, g
