import numpy as np
########################################################################################
# global constants
H = 68.0 # km/s/Mpc
G = 6.67e-20 # km^3/kg s^2
cdmMass = 0.0 # eV
wdmMass = 0.0 # eV
########################################################################################
def integrate():
  return physics

integrate()
# that's it we've solved it
########################################################################################
def gravPotential(x1,x2,dmType): 
  # x1 and x2 are the position vectors of 2 particles; dmType = 0 for WDM, = 1 for CDM
  if dmType == 0:
    m = wdmMass
  else if dmType == 1:
    m = cdmMass
  else 
    print("Invalid dmType parameter given. Must be 0 (WDM) or 1 (CDM).")
    return NaN
  
  r = np.sqrt(((x1[0]-x2[0])**2)+((x1[1]-x2[1])**2)+((x1[2]-x2[2])**2))
  return -1*G*(m**2)/r # this is a scalar -- does it need to be a vector and if so how do we get direction
