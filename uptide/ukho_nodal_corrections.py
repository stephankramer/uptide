# This module implements the nodal corrections based on the 
# "Tidal Harmonic Constants Product Specification" (see: http://www.ukho.gov.uk/AdmiraltyPartners/FGHO/Pages/TidalHarmonics.aspx)
# in particular appendix B

import math
from math import sin,cos,pi

zero_group = ['Z0', 'Sa', 'Ssa', 'pi1', 'P1', 'S1', 'S1', 'S1', 'psi1', 'T2', 'S2', 'R2']

# this covers groups y, a, b c, f, g, j, k, m, and o
group1 = ['M1B', 'M1', 'M1A', 'gamma2', 'alpha2', 'delta2', 'xi2', 'eta2', 'L2']
group2 = ['Mm', 'Mf', 'Mfm', 'O1', 'K1', 'tau1', 'M2', 'K2', 'MSf', 'MSo', 'MSqm', '2SM']
groupf = ['NA2', 'NA2', 'NB2', 'NA2*', 'MA2', 'MB2', 'MA2*'] # same nodal correction as M2
group2 += groupf
groupg = ['M3', 'M5', 'M7'] # odd powers of M (except M1)
group2 += groupg
groupj = ['J1', 'chi1', 'phi1', 'theta1'] # all the same as J1
group2 += groupj
groupm = ['Mqm', 'eps2', 'mu2', 'N2', 'nu2', 'lambda2'] + ['MP1', 'MP1', '2N2'] # also same nodal correction as M2
group2 += groupm
groupo = ['sigma1', 'Q1', 'rho1'] + ['2Q1', 'nuJ1'] # same as O1
group2 += groupo

groupe = ['ups1', 'OO1']

def nodal_corrections(constituents, p, N, pp):
  us = []; fs = []
  for constituent in constituents:
    if constituent in zero_group:
      u,f = 0.,1.
    elif constituent in group1:
      u,f = group1_nodal_correction(constituent, p, N, pp)
    elif constituent in group2:
      u,f = group2_nodal_correction(constituent, p, N, pp)
    else:
      if constituent in ['ups1', 'OO1']: # this is group e
        c = 'KQ1'
      elif constituent=='L2A': # this is group p
        c = '2MN2'
      else:
        c = constituent
      print c

    us.append(u)
    fs.append(f)

  return fs, us

def group1_nodal_correction(constituent, p, N, pp):
  if constituent=='M1B':
    uf = 1j * (2.783*sin(2*p) +  0.558*sin(2*p-N) + 0.184*sin(N)) + \
           1 + 2.783*cos(2*p) +  0.558*cos(2*p-N) + 0.184*cos(N)
  elif constituent=='M1':
    uf = 1j * (sin(p) + 0.2*sin(p-N)) + 2*(cos(p) +  0.2*cos(p-N))
  elif constituent=='M1A':
    uf = 1j * (-0.3593*sin(2*p) - 0.2*sin(N) - 0.066*sin(2*p-N)) + \
	 1  +  0.3593*cos(2*p) +  0.2*cos(N) + 0.066*cos(2*p-N) 
  elif constituent=='gamma2':
    uf = 1j * 0.147*sin(2*(N-p)) + 1. +  0.147*cos(2*(N-p))
  elif constituent=='alpha2':
    uf = 1j * -0.0446*sin(p-pp) + 1. - 0.0446*cos(p-pp)
  elif constituent=='delta2':
    uf = 1j * 0.477*sin(N) + 1. - 0.477*cos(N)
  elif constituent=='xi2' or constituent=='eta2':
    uf = 1j * -0.439*sin(N) + 1. + 0.439*cos(N)
  elif constituent=='L2':
    uf = 1j * (-0.2505*sin(2*p)-0.1102*sin(2*p-N)-0.0156*sin(2*p-2*N)-0.037*sin(N)) + \
             1.-0.2505*cos(2*p)-0.1102*cos(2*p-N)-0.0156*cos(2*p-2*N)-0.037*cos(N)
  else:
    raise Exception("constituent not in group1") 
  return abs(uf), math.atan2(uf.imag, uf.real)

def group2_nodal_correction(constituent, p, N, pp):
  if constituent=='Mm' or constituent=='Mfm': # Mfm is group a
    u = 0.
    f = 1 -  0.1311*cos(N)  +  0.0538*cos(2*p)  +  0.0205*cos(2*p-N)
  elif constituent=='Mf':
    u=-23.7*sin(N)+2.7*sin(2*N)-0.4*sin(3*N)
    f=1.084+0.415*cos(N)+0.039*cos(2*N)
  elif constituent=='O1' or constituent in groupo:
    u=10.80*sin(N)-1.34*sin(2*N)+0.19*sin(3*N)
    f=1.0176+0.1871*cos(N)-0.0147*cos(2*N)
  elif constituent=='K1' or constituent=='tau1': # tau1 is groupk
    u=-8.86*sin(N)+0.68*sin(2*N)-0.07*sin(3*N)
    f=1.0060+0.1150*cos(N)-0.0088*cos(2*N)+0.0006*cos(3*N)
  elif constituent in groupj:
    u=-12.94*sin(N)+1.34*sin(2*N)-0.19*sin(3*N)
    f=1.1029+0.1676*cos(N)-0.0170*cos(2*N)+0.0016*cos(3*N)
  elif constituent=='M2' or constituent in groupf or constituent in groupm:
    u=-2.14*sin(N)
    f=1.0007-0.0373*cos(N)+0.0002*cos(2*N)
  elif constituent=='K2':
    u=-17.74*sin(N)+0.68*sin(2*N)-0.04*sin(3*N)
    f=1.0246+0.2863*cos(N)+0.0083*cos(2*N)-0.0015*cos(3*N)
  elif constituent in groupg:
    s=int(constituent[-1])
    u=-s*1.07*sin(N)
    f=(1.0007-0.0373*cos(N)+0.0002*cos(2*N))**(s/2.)
  elif constituent in ['MSf', 'MSo', 'MSqm']: # group b
    u=2.14*sin(N)
    f=1.0007-0.0373*cos(N)+0.0002*cos(2*N)
  elif constituent=='2SM': # group c
    u=2.*2.14*sin(N)
    f=(1.0007-0.0373*cos(N)+0.0002*cos(2*N))**2.
  else:
    raise Exception("constituent not in group2")
  return u, f

#diurnal_letters=['sigma', 'Q', 'rho', 'O', 'tau', 'K', 'chi', 'pi', 'P', 'psi', 'phi','J', 'ups']
#semidiurnal_letters=['eps', 'mu', 'N', 'nu', 'gamma', 'alpha', 'M', 'delta', 'lambda', 'L', 'T', 'S', 'R', 'xi', 'eta']
#def decompose_constituent(constituent):
