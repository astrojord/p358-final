import numpy as np
from BHA import Body
from analysis import body_outfile

G = 4.5245e-12 #kpc^3/(Msun*My^2)

def NFW_dist(r, rs):
    dens = (r/rs)**(-1.)*(1.+r/rs)**(-2.)
    return dens
def cusp_dist(r, rs):
    dens = r**(-1.)*(r+rs)**(-3.)
    return dens
def NFW_vel(r, rs, M, Rvir):
    v = np.sqrt(3.*G*M*rs**2. / Rvir**3. / (Rvir/rs) / (1.+Rvir/rs)**2. * (np.log(1.+r/rs) / (r/rs) - 1./(1.+r/rs)))
    return v
def cusp_vel(r, rs, M, Rvir):
    v = np.sqrt(G*M/12./rs * (12.*r*(r+rs)**3./rs**4. * np.log((r+rs)/r) - r/(r+rs) * (25.+52.*r/rs+42.*(r/rs)**2.+12.*r/rs)))
    return v

def get_halo(fMH, fDENS, fVELS, rs, delta, npart, Mv, Rv, halonum):
    """
    Generates body object and output file for one initial halo from theoretical density profile
    inputs:
        fMH : function
            function handle for metropolis hastings function that returns array of length npart of r values from sampled from profile
        fDENS : function
            function handle for density
        fVELS : function
            function handle for velocity
        rs : float
            radius scaling factor for initial profile
        delta : float
            MH delta value
        npart : int
            number of particles (in one halo)
        Mv : float
            halo virial mass IN SIMULATION UNITS
        Rv : float
            halo virial radius IN SIMULATION UNITS
        halonum : int
            halo number for the body objects for each particle

    outputs:
        init_halo : list
            list of initial body objects for both halos
    """

    mpp = Mv/npart #mass per particle

    #get r values from methast for each halo
    r0 = 1.
    r = fMH(fDENS,rs,int(npart),r0,delta)

    #uniformly select phi and theta for positions
    phi = np.random.uniform(0,2.*np.pi, len(r))
    theta = np.random.uniform(0,np.pi, len(r))

    #convert to cartesian
    pos = np.array([r*np.sin(theta)*np.cos(phi), r*np.sin(theta)*np.sin(phi), r*np.cos(theta)])

    #calcualte velocity magnitudes
    velmag = fVELS(r, rs, Mv, Rv)

    #velocity direction random in 3D for each particle
    dir = np.random.uniform(size=(npart,3))
    vel = np.zeros((3,npart))
    for i in range(npart):
        dir[i,:] = dir[i,:]/np.linalg.norm(dir[i,:]) #unit vector
        vel[:,i] = velmag[i] * dir[i,:]

    #make a list of body objects to output
    bod = [None] * (npart)
    for i in range(npart):
        bod[i] = Body(pos=pos[:,i],mass=mpp,vel=vel[:,i],acc=np.zeros(3),halonum = halonum)

    #make an output file with the list of body objects
    body_outfile(bod, 'init_halo.txt')

    return bod

