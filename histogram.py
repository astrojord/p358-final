import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

sns.set()

# =============================================================================
# Input: 
# poslist: a list of position vectors 
# Output: a histogram of the radii from the center 
# =============================================================================
def histogram(xarray, yarray, zarray, softrad, n, snapnum): 
    N = len(xarray)
    x = np.sum(xarray)
    y = np.sum(yarray)
    z = np.sum(zarray)
        
    center = np.array((x/N, y/N, z/N))
    
    radii = np.zeros(N)
    for i in range(N):
        xi = xarray[i]
        yi = yarray[i]
        zi = zarray[i]
        xc, yc, zc = center
        
        radii[i] = np.sqrt((xi - xc)**2 + (yi - xc)**2 + (zi - zc)**2) #finds the different radii using the distance formula
        
    binnum = int(np.sqrt(N))
    logradii = np.log10(radii)
    
    hist, edges= np.histogram(logradii,bins = binnum,normed=False)
    #print(hist)
    logx             = 0.5*(edges[0:edges.size-1]+edges[1:edges.size])
    shellvols     = np.zeros(logx.size)
    shellvols[0]  = (4./3.)*np.pi*10**edges[1]**3
    shellvols[1:] = ((4./3.)*np.pi*10**edges[2:]**3)-((4./3.)*np.pi*10**edges[1:edges.size-1]**3)
    hist     = np.divide(hist,shellvols)
    
    x=10**logx
    N = len(hist)
    ymin = min(hist*x**n)#+ 10**(-12) sometimes this is 0 and it won't plot so ya need to add something to make it not zero... not ideal
    ymax = max(hist*x**n)
    print(ymax)
    plt.loglog(x,hist*x**n,'.',[softrad*2.8, softrad*2.8], [ymin, ymax], '-')
    plt.xlabel('log(radius) (kpc)')
    plt.ylabel('log(density*radius**%s)' % str(n))
    plt.title(("Density profile"))
    plt.savefig('snapshot%s.png'% str(snapnum))
    plt.legend(('simulation data', 'softening radius * 2.8'))
    plt.figure
    plt.show()
    
    #ignore below, it doesn't take into account the volume of the shells 
# =============================================================================
#     binnum = int(np.sqrt(N))
#     plt.hist(radii, bins = binnum) 
#     plt.title(('Histogram of radii'))
#     plt.xlabel(('radius'))
#     plt.ylabel(('Counts'))
#     plt.figure
#     plt.show()
# =============================================================================

def extract_data(dataarray):
    N = len(dataarray)
    xarray = np.zeros((N))
    yarray = np.zeros((N))
    zarray = np.zeros((N))
    
    for i in range(N):
        oneparticle = dataarray[i]
        xarray[i] = oneparticle[0]
        yarray[i] = oneparticle[1]
        zarray[i] = oneparticle[2]

    return xarray, yarray, zarray

#for the cuspy (simulation E)
'''
n = 0
snapnum = 119.0 #the .0 means n = 0
dataarraycusp = np.loadtxt(r"snapshot119.txt")
xarray, yarray, zarray = extract_data(dataarraycusp)
softrad = 0.6 #softening radius
print(histogram(xarray, yarray, zarray, softrad, n, snapnum))
'''

'''
n = 5
snapnum = 119.5 #the .5 means n = 5
dataarraycusp = np.loadtxt(r"snapshot119.txt")
xarray, yarray, zarray = extract_data(dataarraycusp)
softrad = 0.6 #softening radius
print(histogram(xarray, yarray, zarray, softrad, n, snapnum))
'''

#for NFW-spindles (simulation A)
'''
n = 1
snapnum = 198
dataarrayNFW = np.loadtxt(r"snapshot%s.txt" % str(snapnum))
xarray, yarray, zarray = extract_data(dataarrayNFW)
softrad = 0.75 #softening radius
print(histogram(xarray, yarray, zarray, softrad, n, snapnum))
'''

#for NFW without spindles (simulation B)
'''
n = 9/4
snapnum = 199.2 #the .2 means n = 9/4
dataarrayNFW2 = np.loadtxt(r"snapshot199.txt")
xarray, yarray, zarray = extract_data(dataarrayNFW2)
softrad = 0.6 #softening radius
print(histogram(xarray, yarray, zarray, softrad, n, snapnum))
'''

'''
n = 0
snapnum = 199.0 #the .0 means n = 0
dataarrayNFW2 = np.loadtxt(r"snapshot199.txt")
xarray, yarray, zarray = extract_data(dataarrayNFW2)
softrad = 0.6 #softening radius
print(histogram(xarray, yarray, zarray, softrad, n, snapnum))
'''

