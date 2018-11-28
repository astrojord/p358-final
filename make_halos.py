import numpy as np

#call methast to get rs and put sph sym other coord; and velocities; make list of body obj; shift to have two halos, id which (thurs)
#assuming 
def get_halos(fMH, npart, profile, sep):
    """
    inputs:
        fMH : function
            function handle for metropolis hastings function that returns array of length npart of r values from sampled from profile
        npart : int
            number of particles
        profile : string
            'infall' or 'NFW'
        sep : float
            separation between centers of halos
    outputs: list
            list of body objects of both halos
    """
    #get rs values from methast for each halo
    #get v^2, then v
    #uniformly select phi and theta for positions
    #put phi hat direction on v
    #convert to cartesian
    #shift one halo sep/2 to the left and the other sep/2 to the right
    #make a list of body objects to output (obj includes position, velocity, and indicator of which halo it is in (initially)
