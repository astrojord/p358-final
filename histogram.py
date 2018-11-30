import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
# Input: 
# poslist: a list of position vectors 
# Output: a histogram of the radii from the center 
# =============================================================================
def histogram(poslist): 
    N = len(poslist)
    x = 0 
    y = 0
    z = 0
    
    for i in range(N):
        posvec = poslist[i] #position vector 
        xi, yi, zi = posvec
        x += xi
        y += yi
        z += zi
        
    center = np.array((x/N, y/N, z/N))
    
    radii = np.zeros(N)
    for i in range(N):
        posvec = poslist[i] #position vector 
        xi, yi, zi = posvec
        xc, yc, zc = center
        
        radii[i] = np.sqrt((xi - xc)**2 + (yi - xc)**2 + (zi - zc)**2) #finds the different radii using the distance formula
        
    binnum = int(np.sqrt(N))
    plt.hist(radii, bins = binnum) 
    plt.title(('Histogram of radii'))
    plt.xlabel(('radius'))
    plt.ylabel(('Counts'))
    plt.figure
    plt.show()
    
