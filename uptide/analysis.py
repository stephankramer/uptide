import numpy
import numpy.linalg

def harmonic_analysis(tide, x, t):
  """Perform tidal harmonic analysis for a given signal x at times t. 
  Returns the amplitudes and phases of the constituents defined in tide 
  (a Tides object), in the order of tide.constituents. The times t are 
  in seconds after the date time set with tide.set_initial_time()."""
  if len(x)!=len(t):
    raise Exception("Length of x and t should be the same")
  N = len(t)
  M = len(tide.omega)

  """We first target for a solution of the form:
    eta = Z0 + \sum_n B_n cos omega_n t + C_n sin omega_n t
        = Z0 + \sum_n Re[ (B_n-i C_n) e^(i omega_n t) ]"""
  if "Z0" in tide.constituents:
    """we solve for the least squares approximation
        y=[B_1,B_2,...B_M,C_1,C_2,...C_M]
      that closest approximates x at all t, where for some j:
      B_j=Z0 and we leave out C_j (as its associated sin() is always zero)"""
    # the indices of the constituents other than Z0
    nonz0 = numpy.array(tide.constituents)!='Z0'
    mat = numpy.hstack([
      numpy.cos(numpy.outer(t,tide.omega)), 
      numpy.sin(numpy.outer(t,tide.omega[nonz0]))
      ])
  else:
    """we solve for the least squares approximation
        y=[Z0,B_1,B_2,...B_M,C_1,C_2,...C_M]
      that closest approximates x at all t"""
    mat = numpy.hstack([
      numpy.ones((N,1)), 
      numpy.cos(numpy.outer(t,tide.omega)), 
      numpy.sin(numpy.outer(t,tide.omega))
      ])

  y = numpy.linalg.lstsq(mat, x)

  if "Z0" in tide.constituents:
    B = y[0][0:M]
    # the C_j associated with Z0 is left zero
    C = numpy.zeros(M)
    C[nonz0] = y[0][M:]
  else:
    Z0 = y[0][0]
    B = y[0][1:M+1]
    C = y[0][M+1:]
  A = B - 1j*C
  """Now with A_n=B_n -i C_n, we have
     eta = Z0 + \sum_n Re[ A_n e^(i omega_n t) ]
         = Z0 + \sum_n Re[ C_n f_n e^(i (-g+phi+u)) e^(i omega_n t)]
         = Z0 + \sum_n a_n f_n cos(omega_n t-g+phi_n+u_n),
  i.o.w. A_n = a_n f_n e^(i (-g+phi+u))"""
  a = numpy.abs(A)/tide.f
  arg = numpy.angle(A)
  g = (tide.phi+tide.u-arg) % (2*numpy.pi)
  return a, g
