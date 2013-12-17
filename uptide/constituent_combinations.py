# module to deal with the mess that is tidal constituent combinations
import numpy

# note that K can be either semi-diurnal or diurnal
diurnal_letters=['SIGMA', 'Q', 'RHO', 'O', 'TAU', 'CHI', 'PI', 'P', 'PSI', 'PHI','J', 'UPS', 'THETA', 'K']
semidiurnal_letters=['EPS', 'MU', 'N', 'NU', 'GAMMA', 'ALPHA', 'M', 'DELTA', 'LAMBDA', 'L', 'T', 'S', 'R', 'XI', 'ETA']

greek_letters=['ALPHA', 'BETA', 'GAMMA', 'DELTA', 'EPS', 'ZETA', 'ETA', 
'THETA', 'LAMBDA', 'MU', 'NU', 'XI', 'PI', 'RHO', 'SIGMA', 'TAU', 'UPS',
'PHI', 'CHI', 'PSI']

# the following combinations occur in the UKHO table, but I can't interpret them:
ukho_numbers_dont_add_up = ["3MS2", "3MS5", "MSP2", "4MS4", "4ML12", "2MNO6", "2MS3", "OO1", "3(SM)N2", "MPS2", "2MS5", "5MSN12", "MS3", "MS1"]
ukho_frequencies_dont_add_up = ["NSK5", "3N2MS12", "M(KS)2", "2NKMS5", "4MSN8", "4M2SN10", "M(SK)2", "5MSN10"]
ukho_unknown_combinations = ukho_numbers_dont_add_up + ukho_frequencies_dont_add_up

class ConstituentCombinationException(Exception):
  pass

# any combination that ends in any of those is long-period and should be covered already
period_letters=['M','A','F','O','N']
# combinations that end in A, B or C
abc_combinations=['M1B', 'M1C', 'M1A', 'L2A', 'L2B']

def decompose_constituents(constituent):
  if constituent in ukho_unknown_combinations:
    # give up straight away on these
    raise ConstituentCombinationException("Don't know how to interpret this constituent combination")
  if constituent in abc_combinations:
    diurnal = int(constituent[1])
    combination = constituent[0]
  elif constituent[0] in ['M','N'] and constituent[1] in ['A','B'] \
      and len(constituent)>2 and constituent[2].isdigit():
    diurnal = int(constituent[2])
    combination = constituent[0] # let's ignore the A or B
  elif constituent=='SNU2':
    diurnal = 0
    combination = 'SNU'
  elif constituent[-1].isdigit():
    if constituent[-2].isdigit():
      diurnal = int(constituent[-2:])
      combination = constituent[:-2]
    else:
      diurnal = int(constituent[-1])
      combination = constituent[:-1]
  elif constituent[-1] in period_letters:
    # don't bother decomposing as all these long-period components are a mess
    # we simply return it as if it were a single constituent
    # (even for combinations like SM, 2SM, SN, 2SMN, etc.)
    return constituent[:-1], [0], [1], 0
  else:
    print constituent
    raise Exception("Unable to interpret constituent combination")

  letters, periods, multiplicities = decompose_combination(combination)

  # try flipping signs, and K1 vs. K2 to ensure sum(periods*multiplicities)==diurnal
  fix_periods_and_signs(letters, periods, multiplicities, diurnal)
  return letters, periods, multiplicities, diurnal

def decompose_combination(combination):
  in_bracket = False
  multiplicity = 1
  letters = []
  multiplicities = []
  periods = []
  while len(combination)>0:
    if combination[0].isdigit():
      multiplicity = int(combination[0])
      combination = combination[1:]
    elif combination[0]=='(':
      in_bracket = True
      combination = combination[1:]
    elif combination[0]==')':
      in_bracket = False
      multiplicity = 1
      combination = combination[1:]
    else:
      for x in greek_letters:
        if combination.startswith(x):
          letter = x
          combination = combination[len(x):]
          break
      else: # not starting with any greek letter
        if not combination[0].isalpha():
          print constituent, combination
          raise Exception("Unexpected character in constituent combination")
        letter = combination[0]
        combination = combination[1:]
      # add the greek or latin letter:
      letters += [letter]
      multiplicities += [multiplicity]
      if letter in diurnal_letters:
        periods.append(1)
      elif letter in semidiurnal_letters:
        periods.append(2)
      elif letter=='Z':
        periods.append(0)
      else:
        print constituent, letter
        raise Exception("Unexpected letter in constituent combination")
      if not in_bracket:
        multiplicity = 1

  return letters, numpy.array(periods), numpy.array(multiplicities)

def fix_periods_and_signs(letters, periods, multiplicities, diurnal):
  """Try to flip the signs of multiplicities and changing between K2 and K1
  to make the sum(periods*multiplicities)==diurnal. 
  Raises ConstituentCombinationException if failed."""
  if len(letters)==1:
    if periods[0]==0 and diurnal==0:
      pass
    elif periods[0]==1:
      multiplicities[0] = diurnal
    elif periods[0]==2:
      if diurnal%2==0:
        multiplicities[0]=diurnal/2
      else:
        periods[0]=1
        multiplicities[0]=diurnal
    else:
      raise ConstituentCombinationException("Unable to interpret constituent combination")
    return

  # try fixing the sum by flipping signs in multiplicities first
  try:
    fix_signs(letters, periods, multiplicities, diurnal)
  except ConstituentCombinationException:
    # then, try flipping the period of K1 to K2
    # start from beginning with signs all positive
    multiplicities[:] = abs(multiplicities)
    for k in range(len(letters)):
      if letters[k]=='K':
        periods[k] = 2
        fix_signs(letters, periods, multiplicities, diurnal)
        return
    else: # no Ks
      raise

def fix_signs(letters, periods, multiplicities, diurnal):
  """Try to flip the signs of multiplicities 
  to make the sum(periods*multiplicities)==diurnal. 
  Raises ConstituentCombinationException if failed."""
  if sum(periods*multiplicities)==diurnal:
    return

  for k in range(len(letters)-1,0,-1):
    multiplicities[k] *= -1
    if sum(periods*multiplicities)==diurnal:
      return

  raise ConstituentCombinationException("No combination of signs adds up for this constituent combination")
