import numpy as np
from BHA import Body

G = 4.5245e-12 #kpc^3/(Msun*My^2)

def NFW_dist(r, rs):
    dens = (r/rs)**(-1.)*(1+r/rs)**(-2.)
    return dens
def cusp_dist(r, rs):
    dens = r**(-1.)*(r+rs)**(-3.)
    return dens
def NFW_vel(r, rs, M, Rvir):
    v = np.sqrt(3.*G*M*rs**2. / Rvir**3. / (Rvir/rs) / (1.+Rvir/rs)**2. * (np.log(1.+r/rs) / (r/rs) - 1./(1.+r/rs)))
    return v
def cusp_vel(r, rs, M, Rvir):
    v = np.sqrt(G*M/12./a * (12.*r*(r+a)**3./a**4. * np.log((r+a)/r) - r/(r+a) * (25.+52.*r/a+42.*(r/a)**2.+12.*r/a)))
    return v

def get_halos(fMH, fDENS, fVELS, rs, delta, npart, M_typ, Rv_typ, bigrat, sep, axis):
    """
    inputs:
        fMH : function
            function handle for metropolis hastings function that returns array of length npart of r values from sampled from profile
        fDENS : function
            function handle for density
        fVELS : function
            function handle for velocity squared
        rs : float
            radius scaling factor for initial profile
        delta : float
            MH delta value
        npart : int
            number of particles (or close to what it will be, after ratios)
        M_typ : float
            typical dark matter galaxy halo mass IN SIMULATION UNITS
        Rv_typ : float
            typical dark matter galaxy halo mass IN SIMULATION UNITS
        bigrat : float
            must be less than 1, ratio of total particles in larger halo (ratio of mass to total mass where Mtotal = 2*Mtyp)
        sep : float
            separation between centers of halos
        axis : bool
            true for separation on axis in plane of rotation, false for separation on axis perpendicular to plane of rotation

    outputs:
        init_halos : list
            list of initial body objects for both halos
    """
    #use bigrat to get M1, M2, npart1, npart2, Rv1, Rv2
    mpp = 2.*M_typ/npart #mass per particle
    M1 = bigrat*2.*M_typ
    M2 = (1.-bigrat)*2.*M_typ
    npart1 = int(bigrat*npart)
    npart2 = int((1-bigrat)*npart)
    Rvir1 = (M1/M_typ)**(1./3.) * Rv_typ  #keeping total density constant
    Rvir2 = (M2/M_typ)**(1./3.) * Rv_typ

    #get r values from methast for each halo
    r10 = 1.
    r20 = 1.
    r1 = fMH(fDENS,rs,int(npart1),r10,delta)
    r2 = fMH(fDENS,rs,int(npart2),r20,delta)

    #calcualte velocities
    vel1 = fVELS(r1, rs, M1, Rvir1)
    vel2 = fVELS(r2, rs, M2, Rvir2)

    #uniformly select phi and theta for positions
    phi1 = np.random.uniform(0,2.*np.pi, len(r1))
    theta1 = np.random.uniform(0,np.pi, len(r1))
    phi2 = np.random.uniform(0,2.*np.pi, len(r2))
    theta2 = np.random.uniform(0,np.pi, len(r2))

    #convert to cartesian
    pos1 = np.array([r1*np.sin(theta1)*np.cos(phi1), r1*np.sin(theta1)*np.sin(phi1), r1*np.cos(theta1)])
    pos2 = np.array([r2*np.sin(theta2)*np.cos(phi2), r2*np.sin(theta2)*np.sin(phi2), r2*np.cos(theta2)])

    vel1 = vel1*np.array([-1.*np.sin(phi1), np.cos(phi1), np.zeros(npart1)])
    vel2 = vel2*np.array([-1.*np.sin(phi2), np.cos(phi2), np.zeros(npart2)])

    #shift one halo sep/2 to the left and the other sep/2 to the right
    if axis == True:
        pos1[0,:] += sep/2.
        pos2[0,:] -= sep/2.
    elif axis == False:
        pos1[2,:] += sep/2.
        pos2[2,:] -= sep/2.

    #make a list of body objects to output
    bod = [None] * (npart1 + npart2)
    for i in range(npart1):
        bod[i] = Body(pos=pos1[:,i],mass=mpp,vel=vel1[:,i],acc=np.zeros(3),halonum = 1)
    for i in range(npart2):
        bod[npart1 + i] = Body(pos=pos2[:,i],mass=mpp,vel=vel2[:,i],acc=np.zeros(3),halonum = 2)
    return bod

