import numpy as np
from BHA import Body

def gen_halos(filepath,sep):

    # x y z m vx vy vz ax ay az haloNum
    # 0 1 2 3 4  5  6  7  8  9  10

    data = np.loadtxt(filepath) # insert file path here xd
    data[:,0] = data[:,0] - sep/2
    data[:,4] = data[:,4] + 3 # kpc/Myr

    data2 = data.copy()
    data2[:,10] = 2 # change halonum to 2
    data2[:,0] = data2[:,0] + sep/2
    data2[:,4] = data2[:,4] - 3 # kpc/Myr

    bodies = []
    nP = data.shape[0]
    for i in range(nP):
        bodies.append(Body(pos = data[i,0:3], mass = data[i,3],vel = data[i,4:7], acc = data[i,7:10], halonum = 1))
        bodies.append(Body(pos = data2[i,0:3], mass = data2[i,3],vel = data2[i,4:7], acc = data2[i,7:10], halonum = 2))


    print(bodies)
    return bodies
