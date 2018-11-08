import numpy as np
from BHA import Body
########################################################################################
# global constants
H = 68.0  # km/s/Mpc
G = 6.67e-20  # km^3/kg s^2

omegaM = 0  # matter in total (omegaB + omegaDM)
omegaB = 0  # baryonic matter
omegaDM = 0  # dark matter
omegaDE = 0  # dark energy

cdmMass = 0.0  # eV
wdmMass = 0.0  # eV
########################################################################################
def get_dxdt(bod, tau, ptcl_tree):
    """
    calculates total acceleration from equation of motion (yall can do this for
    one body at a time and put the loop in integrate, or just all the bodies)

    inputs
    -----------------
    bod : Body
      one individual particle with ndims position
    neighbor_tree : Tree
      the tree structure containing neighbor CM and M data

    outputs
    -----------------
    accel : ndarray
      3x1 array of dx_i/dtau
    """
    assert type(bod) == Body, "bod input must be a Body object"
    # get neighbor list
    neighbor_list = ptcl_tree.neighbors(bod)
    
    # calculate the negative of the gravitational potential
    negGradPotential = np.zeros(3)
    for neigh in neighbor_list:
        posit = neigh[0]
        mass = neigh[1]
        dvect = bod.pos
        r = posit - dvect
        for i in range(0,2):
            negGradPotential[i] = -1*(G*mass*r[i])/((np.dot(r,r))**(3/2))

    # calculate the H(tau)v(tau) term
    expansion = np.zeros(3)
    a = a(tau,0)
    # how to get v(t)!

    accel = negGradPotential - expansion
    return accel

def a(time, mode):
    """
    get a (cosmological scale factor) in terms of either conformal or
    cosmological time (tau or t)

    equations from Sazhin et. al. 2011

    inputs
    ----------------
    time : float
      value of tau or t, depending on mode
    mode : int
      0 or 1; 0 corresponds to conformal time tau, 1 corresponds to
              cosmological time t
    outputs
    ----------------
    a : float
    """
    if mode == 0:  # conformal
        return (omegaM / 4) * (H0 * time) ** 2
    else if mode == 1:  # cosmological
        return (9 * omegaM / 4) ** (1 / 3) * (H0 * time) ** (2 / 3)
    else
        print("invalid mode input -- must be 0 or 1")
        return Nan

########################################################################################
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
# bo = Body((3,4,5), 60)
# print(bo)
########################################################################################
