# This module implements the nodal corrections based on the 
# "Tidal Harmonic Constants Product Specification" (see: http://www.ukho.gov.uk/AdmiraltyPartners/FGHO/Pages/TidalHarmonics.aspx)
# in particular appendix B

import math
from math import sin,cos,pi
import constituent_combinations

class NodalCorrectionException(Exception):
  pass

# group with no nodal corrections
zero_group = ['Z0', 'SA',  'SSA', 'PI1', 'P1', 'S1', 'S1', 'S1', 'PSI1', 'T2', 'S2', 'R2']
unknown_group = ['STA', 'MSTM']

group1 = ['M1B', 'M1', 'M1A', 'GAMMA2', 'ALPHA2', 'DELTA2', 'XI2', 'ETA2', 'L2']
group2 = ['MM', 'MF', 'MFM', 'O1', 'K1', 'TAU1', 'M2', 'K2', 'MSF', 'SMF', 'MS0', 'MSQM', '2SM', 'LA2']
groupe = ['UPS1', 'OO1'] # same nodal correction as KQ1
group2 += groupe
groupf = ['NA2', 'NA2', 'NB2', 'NA2*', 'MA2', 'MB2', 'MA2*'] # same nodal correction as M2
group2 += groupf
groupg = ['M3', 'M5', 'M7'] # odd powers of M (except M1)
group2 += groupg
groupj = ['J1', 'CHI1', 'PHI1', 'THETA1'] # all the same as J1
group2 += groupj
groupm = ['MQM', 'EPS2', 'MU2', 'N2', 'NU2', 'LAMBDA2'] + ['MP1', 'MP1', '2N2'] + ['MSm', 'SM'] # also same nodal correction as M2
group2 += groupm
groupo = ['SIGMA1', 'Q1', 'RHO1'] + ['2Q1', 'NUJ1'] # same as O1
group2 += groupo
# this group I've made up myself:
# I presume KO0=K1-O1
# MK seems to be the same as KO0
# SN is S2-N2 and SM is S2-M2, but since nc of N2==nc of M2 and nc of S==0, we get nc of SN=-nc of M2
# 2SMN=SM-SN=2S2-M2-N2 ?
groupq = ['KO0', 'MK0', 'SN', 'SM']
group2 +=groupq

groupe = ['UPS1', 'OO1']

def nodal_corrections(constituents, p, N, pp):
  fs = []; us = []
  for constituent in constituents:
    try:
      f, u = nodal_corrections_single_constituent(constituent, p, N, pp)
    except NodalCorrectionException:
      f = 1.0; u =0.0
      letters, periods, multiplicities, diurnal = constituent_combinations.decompose_constituents(constituent)
      for l,p,m in zip(letters, periods, multiplicities):
        fi, ui = nodal_corrections_single_constituent(l+str(p), p, N, pp)
        f *= fi**abs(m)
        u += m*ui
    fs.append(f)
    us.append(u)

  return fs, us


def nodal_corrections_single_constituent(constituent, p, N, pp):
  constituent = constituent.upper()
  if constituent[-1]=='O':
    constituent[-1]='0'
  if constituent in zero_group:
    u,f = 0.,1.
  elif constituent in group1:
    u,f = group1_nodal_correction(constituent, p, N, pp)
  elif constituent in group2:
    u,f = group2_nodal_correction(constituent, p, N, pp)
  else:
    print constituent
    raise NodalCorrectionException("This is not a single constituent")

  return f, u

def group1_nodal_correction(constituent, p, N, pp):
  if constituent=='M1B':
    uf = 1j * (2.783*sin(2*p) +  0.558*sin(2*p-N) + 0.184*sin(N)) + \
           1 + 2.783*cos(2*p) +  0.558*cos(2*p-N) + 0.184*cos(N)
  elif constituent=='M1':
    uf = 1j * (sin(p) + 0.2*sin(p-N)) + 2*(cos(p) +  0.2*cos(p-N))
  elif constituent=='M1A':
    uf = 1j * (-0.3593*sin(2*p) - 0.2*sin(N) - 0.066*sin(2*p-N)) + \
	 1  +  0.3593*cos(2*p) +  0.2*cos(N) + 0.066*cos(2*p-N) 
  elif constituent=='GAMMA2':
    uf = 1j * 0.147*sin(2*(N-p)) + 1. +  0.147*cos(2*(N-p))
  elif constituent=='ALPHA2':
    uf = 1j * -0.0446*sin(p-pp) + 1. - 0.0446*cos(p-pp)
  elif constituent=='DELTA2':
    uf = 1j * 0.477*sin(N) + 1. - 0.477*cos(N)
  elif constituent=='XI2' or constituent=='ETA2':
    uf = 1j * -0.439*sin(N) + 1. + 0.439*cos(N)
  elif constituent=='L2':
    uf = 1j * (-0.2505*sin(2*p)-0.1102*sin(2*p-N)-0.0156*sin(2*p-2*N)-0.037*sin(N)) + \
             1.-0.2505*cos(2*p)-0.1102*cos(2*p-N)-0.0156*cos(2*p-2*N)-0.037*cos(N)
  else:
    raise NodalCorrectionException("constituent not in group1") 
  return abs(uf), math.atan2(uf.imag, uf.real)

def group2_nodal_correction(constituent, p, N, pp):
  if constituent=='MM' or constituent=='MFM': # Mfm is group a
    u = 0.
    f = 1 -  0.1311*cos(N)  +  0.0538*cos(2*p)  +  0.0205*cos(2*p-N)
  elif constituent=='MF':
    u=-23.7*sin(N)+2.7*sin(2*N)-0.4*sin(3*N)
    f=1.084+0.415*cos(N)+0.039*cos(2*N)
  elif constituent=='O1' or constituent in groupo:
    u=10.80*sin(N)-1.34*sin(2*N)+0.19*sin(3*N)
    f=1.0176+0.1871*cos(N)-0.0147*cos(2*N)
  elif constituent=='K1' or constituent=='TAU1': # tau1 is groupk
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
  elif constituent=='SM' or constituent=='SN': # basically -M2 (see above groupq)
    u=2.14*sin(N)
    f=1.0007-0.0373*cos(N)+0.0002*cos(2*N)
  elif constituent=='KO0' or constituent=='MK0': # should this be K1-O1?
    # this is a copy from K1:
    u=-8.86*sin(N)+0.68*sin(2*N)-0.07*sin(3*N)
    f=1.0060+0.1150*cos(N)-0.0088*cos(2*N)+0.0006*cos(3*N)
    # subtract O1
    u-=10.80*sin(N)-1.34*sin(2*N)+0.19*sin(3*N)
    f*=1.0176+0.1871*cos(N)-0.0147*cos(2*N)
  elif constituent in groupe: # same as KQ1=K2-Q1=K2-O1
    # copy from K2
    u=-17.74*sin(N)+0.68*sin(2*N)-0.04*sin(3*N)
    f=1.0246+0.2863*cos(N)+0.0083*cos(2*N)-0.0015*cos(3*N)
    # subtract O1
    u-=10.80*sin(N)-1.34*sin(2*N)+0.19*sin(3*N)
    f*=1.0176+0.1871*cos(N)-0.0147*cos(2*N)
  elif constituent=='L2A': # groupp: same as 2MN2=2*M2-N2
    # u_2MN2=2*u_M2-u_N2=u_M2
    u=-2.14*sin(N)
    # f_2MN2=f_M2**2.*f_N2=(f_M2)**3.
    f=(1.0007-0.0373*cos(N)+0.0002*cos(2*N))**3.
  else:
    print constituent
    raise NodalCorrectionException("constituent not in group2")
  return u, f
