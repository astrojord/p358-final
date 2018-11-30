import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation

data = np.loadtxt('file.txt') # assuming just a list or something of data values for each particle
nParticles = 1 # ?
timeSteps = 1000 # ?

# this animates the position of 1 particle given a 3D array of the x,y,z coord (dimensions timeSteps,timeSteps,timeSteps)
pointData = np.array([[1,0,4],[2,2,7],[3,5,1]]) # in reality this will probably be a splice of the data array

fig = plt.figure()
ax = p3.Axes3D(fig)

def updatePoint(num, pointData, lines):
    for line, data in range(0,timeSteps):
        line.set_data(pointData[:,:,num]) # set x and y
        line.set_3d_properties(data[num,num,:]) # set z
    return lines
    
ax.set_xlim3d([-15.0, 15.0])
ax.set_xlabel('x (kpc)')

ax.set_ylim3d([-15.0, 15.0])
ax.set_ylabel('y (kpc)')

ax.set_zlim3d([-15.0, 15.0])
ax.set_zlabel('z (kpc)')

ax.set_title("Smashing halos: animation")

anim = animation.FuncAnimation(fig, updatePoint, timeSteps)
fig.show()

# to do: extend to multiple particles without running 1 trillion different functions for 1 trillion time steps
#        make sure the update function actually works properly in terms of updating x/y and z separately
