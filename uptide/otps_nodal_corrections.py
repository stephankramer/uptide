import numpy
import math

pi = math.pi
deg2rad = pi/180

nodal_correction_f0 = {
  'Mf': 1.043, 
  'O1': 1.009,
  'Q1': 1.009,
  'K1': 1.006,
  'K2': 1.024}
nodal_correction_f1 = {
  'Mm': -0.130,
  'Mf': +0.414,
  'O1': +0.187,
  'Q1': +0.187,
  'K1': +0.115,
  'M2': -0.037,
  'N2': -0.037,
  'K2': +0.286}
nodal_correction_u1 = {
  'Mf': -0.41364303,
  'O1':  0.18849556,
  'Q1':  0.18849556,
  'K1': -0.1553343,
  'M2': -0.03665191,
  'N2': -0.03665191,
  'K2': -0.30892328}
nodal_correction_f2={}
# nodal corrections for M4, MN4, MS4
for comp in ('M2','N2','S2'):
  if comp[0]=='M':
    name='M4'
  else:
    name='M'+comp[0]+'4'
  nodal_correction_f0[name] = nodal_correction_f0.get('M2',1.0) * nodal_correction_f0.get(comp, 1.0)
  nodal_correction_f1[name] = (nodal_correction_f0.get('M2',1.0) * nodal_correction_f1.get(comp, 0.0) +
            nodal_correction_f1.get('M2',0.0) * nodal_correction_f0.get(comp, 1.0))
  nodal_correction_f2[name] = nodal_correction_f1.get('M2',0.0) * nodal_correction_f1.get(comp, 0.0)
  nodal_correction_u1[name] = nodal_correction_u1.get('M2',0.0) + nodal_correction_u1.get(comp, 0.0)
# nodal corrections that are the same as M2 and N2 (see Pugh table 4.3):
for comp in ('2N2', 'MU2', 'NU2', 'T2'):
  nodal_correction_f0[comp] = nodal_correction_f0.get('M2', 1.0)
  nodal_correction_f1[comp] = nodal_correction_f1.get('M2', 0.0)
  nodal_correction_f2[comp] = nodal_correction_f2.get('M2', 0.0)
  nodal_correction_u1[comp] = nodal_correction_u1.get('M2', 0.0)

def nodal_corrections(constituents, N, pp):
  # compute the 18.6 year variations (pp is currently not used)
  # the numbers come from Kowalik and Luick, table 1.6
  f=[]; u=[]
  cosN = math.cos(N*deg2rad)
  cosNsq = cosN**2
  sinN = math.sin(N*deg2rad)
  for constituent in constituents:
    # amplitude corrections:
    f0 = nodal_correction_f0.get(constituent, 1.0)
    f1 = nodal_correction_f1.get(constituent, 0.0)
    f2 = nodal_correction_f2.get(constituent, 0.0)
    f.append(f0 + f1*cosN + f2*cosNsq)
    # phase corrections:
    u.append(nodal_correction_u1.get(constituent, 0.0)*sinN)

  return numpy.array(f), numpy.array(u)
