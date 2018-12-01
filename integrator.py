import numpy as np
from BHA import Body
from BHA import Node
from analysis import body_outfile
########################################################################################
# global constants
H = 68.0  # km/s/Mpc
G = 4.5245e-12 #kpc^3/(Msun*My^2)

omegaM = 0  # matter in total (omegaB + omegaDM)
omegaB = 0  # baryonic matter
omegaDM = 0  # dark matter
omegaDE = 0  # dark energy

cdmMass = 0.0  # eV
wdmMass = 0.0  # eV
########################################################################################
def grav_accelerate(bod, ptcl_tree,n,snap):
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
    snap : int
        1 if a snapshot is being taken, 0 otherwise
    """
    # get neighbor list
    neighbor_list = ptcl_tree.neighbors(bod)
    dvect = bod.pos
    mbod  = bod.mass
    eps = 2**2  #grav softening
    #l_soft2 = 1**2 #square of distance below which softening is used

    accel = np.zeros(len(dvect))
    Ug = 0 #gravitational potential, only actually calculated for snapshots

    for neigh in neighbor_list:

        posit = neigh[0]
        mass = neigh[1]
        d = (posit - dvect)
        dmag2 = np.dot(d,d)
        if dmag2 > 1e-5: #if distance is close to 0 within machine precision they are probably both the same particle
            accel += G * mass / ((dmag2+eps)**1.5) * d

            if snap == 1:
                Ug += -0.5*G*mass*mbod / (np.sqrt(dmag2+eps)) #factor of 1/2 from double counting particle pairs

    return accel, Ug
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

def leapfrog(bods,h,n,l,snap):

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
    snap : int
        1 if a snapshot is being taken, 0 else
    """

    ptcl_tree = Node(pos = np.zeros(n), length = l) #make particle tree

    for bod in bods: #fitting body to tree
        ptcl_tree.fit(bod)

    ptcl_tree.calculate_coms() #calculate coms

    U_total = 0 #total gravitational potential, only calculated for snapshots

    for bod in bods: #update pos/vel/acc of each body with leapfrog equations
        v1 = bod.vel #1 n, 2 is n+1/2, 3 is n+1
        x1 = bod.pos
        F1 = bod.acc

        v2 = v1 + 0.5*h*F1
        x3 = x1 + h*v2
        bod.update_pos(x3) #update pos in bods list
        F3, Ug, = grav_accelerate(bod,ptcl_tree,3,snap)

        v3 = v2 + 0.5*h*F3

        if snap == 1:
            U_total += Ug

        bod.update_vel(v3) #update body attributes in bods list
        bod.update_acc(F3)


    return bods, U_total

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
    nstep = int((tf-ti) // h)
    snapStep = nSave #number of steps between snapshots
    snaps  = []
    U_list = []
    T_list = []
    snapn = 0

    for i in range(nstep):
        print(i)
        if i%snapStep == 0:
            bods, U_total = leapfrog(bods,h,N,l,1)
            snaps.append(np.copy(bods))
            U_list.append(U_total)
            body_outfile(bods,"snapshot%03d.txt"%snapn)
            snapn += 1

            T_total = 0

            for bod in bods:
                v2 = np.dot(bod.vel,bod.vel)
                m = bod.mass
                T = 0.5*m*v2
                T_total += T
            T_list.append(T_total)

        else:
            bods, U_total = leapfrog(bods,h,N,l,0)

    return snaps, U_list, T_list

# use this to see what a body object is
# bo = Body((3,4,5), 60)
# print(bo)
########################################################################################
