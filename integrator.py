import numpy as np
from BHA import Body
from BHA import Node
########################################################################################
# global constants
H = 68.0  # km/s/Mpc
G = 0.8#6.67e-20  # km^3/kg s^2

omegaM = 0  # matter in total (omegaB + omegaDM)
omegaB = 0  # baryonic matter
omegaDM = 0  # dark matter
omegaDE = 0  # dark energy

cdmMass = 0.0  # eV
wdmMass = 0.0  # eV
########################################################################################
def grav_accelerate(bod, ptcl_tree,n):
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
    # get neighbor list
    neighbor_list = ptcl_tree.neighbors(bod)
    dvect = bod.pos
    eps = 0.01 #grav softening

    accel = np.zeros(len(dvect))
    for neigh in neighbor_list:
        posit = neigh[0]
        mass = neigh[1]
        d = posit - dvect
        dmag2 = np.dot(d,d)
        accel += G * mass / ((dmag2+eps)**1.5) * d

    return accel
'''
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

    if mode == 0:  # conformal
        return (omegaM / 4) * (H0 * time) ** 2
    elif mode == 1:  # cosmological
        return (9 * omegaM / 4) ** (1 / 3) * (H0 * time) ** (2 / 3)
    else
        print("invalid mode input -- must be 0 or 1")
        return Nan

'''

########################################################################################


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
        eps = 0.01 #temporary softening term (eventually should only be put in for close encounters)
        for i in range(0,2):
            negGradPotential[i] = -1*(G*mass*r[i])/((np.dot(r,r)+eps)**(3/2))
    # calculate the H(tau)v(tau) term
    #expansion = np.zeros(3)
    #a = a(tau,0)
    # how to get v(t)!
    accel = negGradPotential #- expansion
    return accel
    '''
    outputs
    ----------------
    var : type
      description
    """
    return
    """
    '''

def leapfrog(bods,h,n,l):

    """
    implement one time step of leapfrog to get new positions and
    velocities for all particles
    inputs
    ----------------
    bods : list
      list of body objects, one object for each particle
    h : float
      size of timestep
    n : int
      number of dimensions
    l : flt
      side length of box (units?)
    """

    ptcl_tree = Node(pos = np.zeros(n), length = l) #make particle tree

    for bod in bods: #fitting body to tree
        ptcl_tree.fit(bod)

    ptcl_tree.calculate_coms() #calculate coms

    for bod in bods: #update pos/vel/acc of each body with leapfrog equations
        v1 = bod.vel #1 n, 2 is n+1/2, 3 is n+1
        x1 = bod.pos
        F1 = bod.acc

        v2 = v1 + 0.5*h*F1
        x3 = x1 + h*v2
        bod.update_pos(x3) #update pos in bods list
        F3 = grav_accelerate(bod,1,ptcl_tree)
        v3 = v2 + 0.5*h*F3

        bod.update_vel(v3) #update body attributes in bods list
        bod.update_acc(F3)

    return bods

def integrate(bods_init,ti,tf,h,N,l,nSave):
    '''

    bods_init : list
        list of body objects corresponding to initial particle conditions
    ti : float
        simulation start time
    tf : float
        simulation end time
    h : float
        size of timestep
    N : int
        number of physical dimensions
    l : float
        physical size of box
    nSave : int
        number of snapshots
    '''

    bods = bods_init
    t = ti
    nstep = (tf-ti) // h
    snapStep = nstep // nSave #number of steps between snapshots
    snaps = []

    for i in range(nstep):
        print(i)
        bods = leapfrog(bods,h,N,l)
        if i%snapStep == 0:
            snaps.append(bods)

    return snaps

# use this to see what a body object is
# bo = Body((3,4,5), 60)
# print(bo)
########################################################################################
