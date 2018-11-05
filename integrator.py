import numpy as np
from BHA import Body

########################################################################################
# global constants
H = 68.0 # km/s/Mpc
G = 6.67e-20 # km^3/kg s^2
cdmMass = 0.0 # eV
wdmMass = 0.0 # eV
########################################################################################

def grav_accelerate(bod, ptcl_tree):
    """
    calculates gravitational acceleration on one particle

    inputs
    -----------------
    bod : Body
      one individual particle with ndims position
    neighbor_tree : Tree
      the tree structure containing neighbor CM and M data

    outputs:
    -----------------
    accel : ndarray
      ndims x 1 array of accelerations
    """
    assert type(bod) == Body, "bod input must be a Body object"
    # get neighbor list
    neighbor_list = ptcl_tree.neighbors(bod)

    accel = 0
    for neigh in neighbor_list:
        posit = neigh[0]
        mass = neigh[1]
        dvect = bod.pos
        d = posit - dvect
        dmag2 = np.dot(d,d)
        accel += G * mass / dmag2**1.5 * d

    return accel

def get_dxdt():
    """
    calculates total acceleration from equation of motion (yall can do this for
    one body at a time and put the loop in integrate, or just all the bodies)

    inputs
    -----------------
    var : type
      description
    var : type
      description

    outputs
    -----------------
    var : type
      description
    """
    return


def leapfrog():
    """
    implement the time step equations from leapfrog to get new positions and
    velocities for all particles

    inputs
    ----------------
    var : type
      description

    outputs
    ----------------
    var : type
      description
    """
    return


# use this to see what a body object is
#bo = Body((3,4,5), 60)
#print(bo)












########################################################################################










"""
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
"""
