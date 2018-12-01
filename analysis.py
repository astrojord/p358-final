import numpy as np

def body_outfile(bodylist, saveas):
    """
    saves an output file for data in a list of 3D body objects
    not generalized for n-dims, only 3-space

    inputs:
    body : list
        list of body objects you want to save
    N : int
        number of dimensions (probably 3 tbh)
    saveas : string
        what you want to call the file (and any additonal directory path)

    output:
    header : string
        the name of column headers in order (this will not be in the file,
        in an effor to facilitate easy reading with np.loadtxt(saveas))

    """
    bods = len(bodylist)
    position = np.zeros((bods,3))
    velocity = np.zeros((bods,3))
    accel = np.zeros((bods,3))
    mass = np.zeros(bods)
    halonum = np.zeros(bods)

    for bodn in range(bods):
        position[bodn,:] = bodylist[bodn].pos
        velocity[bodn,:] = bodylist[bodn].vel
        accel[bodn,:] = bodylist[bodn].acc
        mass[bodn] = bodylist[bodn].mass
        halonum[bodn] = bodylist[bodn].halonum

    arrayay = np.array([position[:,0], position[:,1], position[:,2], mass[:],
                        velocity[:,0], velocity[:,1], velocity[:,2],
                        accel[:,0], accel[:,1], accel[:,2], halonum[:]])
    arrayay = np.transpose(arrayay)
    header = ['x', 'y', 'z', 'm', 'vx', 'vy', 'vz', 'ax', 'ay', 'az', 'halonum']

    np.savetxt(saveas,arrayay)

    return header
