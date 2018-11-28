from integrator import integrate
import numpy as np
from BHA import Body
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

N = 3 #number of physical dimensions (plotting set for 3)

nP = 1000 #number of particles
bodz = []
for i in range(nP): #generate some random particles
    pos = (0.5-np.random.rand(N))*10**4
    print(pos)
    vel = (0.5-np.random.rand(N))*0
    bod = Body(pos=pos,mass=1e9,vel=vel,acc=np.zeros(N))
    print(bod)
    bodz.append(bod)
print(bodz)

ti = 0 #initial time
tf = 9000 #final time
h  = 20 #timestep
l = 1e7 #box size, real fuckin big rn because things are flying all over the place
nSave = 10 #number of steps between each save

bodlist = integrate(bodz,ti,tf,h,N,l,nSave)
print(bodz)
pos1 = np.zeros((nP,N))

for i in range(nP):
    for j in range(N):
        pos = bodz[i].pos
        print(bodz[i])
        pos1[i,j] = pos[j]

pos1x = pos1[:,0]
pos1y = pos1[:,1]
pos1z = pos1[:,2]

fig = plt.figure()
ax = fig.add_subplot(111,projection='3d')
ax.plot(pos1x,pos1y,pos1z,'b.')
plt.show()

np.savetxt('dmpos.txt',pos1)

#plt.plot(pos1x,pos1y,'b.')
#plt.show()
