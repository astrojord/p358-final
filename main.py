import numpy as np
import make_halos as init
import methast as mh
import argparse
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from integrator import integrate
import time

# parse arguments
parser = argparse.ArgumentParser()

parser.add_argument('profile', type=str, help='string of the profile type: "NFW" or "cusp"')
parser.add_argument('npart', type=int, help='target total number of particles; note: this might change by +/-1')
parser.add_argument('M', type=float, help='typical galactic dark matter halo mass')
parser.add_argument('R', type=float, help='typical galactic dark matter halo virial radius')
parser.add_argument('rs', type=float, help='scale radius for dark matter halo profiles')
parser.add_argument('delta', type=float, help='delta value to use in metropolis hastings')
parser.add_argument('sep', type=float, help='separation between halos')
parser.add_argument('axis', type=bool, help='true for sep on axis in plane of rotation, false for sep on axis perpendicular to rotation')

args = parser.parse_args()

color = ['b', 'g', 'b'] #for plotting different halos, someone indexed from 1 so in case we ever forget that there are three colors here

#initialize the problem
if args.profile == 'NFW':
    fDENS = init.NFW_dist
    fVELS = init.NFW_vel
elif args.profile == 'cusp':
    fDENS = init.cusp_dist
    fVELS = init.cusp_vel

# get initial halo
halonum = 1
halos_init = init.get_halo(mh.methast, fDENS, fVELS, args.rs, args.delta, args.npart, args.M, args.R, halonum)
bodz = halos_init

#plot and check
pos1 = np.zeros((len(bodz),3))
halo = np.zeros(len(bodz))
for i in range(len(bodz)):
    for j in range(3):
        pos1[i,j] = bodz[i].pos[j]

    halo[i] = bodz[i].halonum

fig = plt.figure()
ax = fig.add_subplot(111,projection='3d')
halo1_ind = np.where(halo == 1)[0]
halo2_ind = np.where(halo == 2)[0]
print(len(halo1_ind))
print(len(halo2_ind))
print(pos1[halo1_ind])
ax.plot(pos1[halo1_ind,0],pos1[halo1_ind,1],pos1[halo1_ind,2], 'r.', markersize=2)
ax.plot(pos1[halo2_ind,0],pos1[halo2_ind,1],pos1[halo2_ind,2], 'b.', markersize=2)
lims = (ax.get_xlim(),ax.get_ylim(),ax.get_zlim())
lim_diffs = (lims[0][1]-lims[0][0],lims[1][1]-lims[1][0],lims[2][1]-lims[2][0])
maxlimInd = np.where(lim_diffs == np.max(lim_diffs))
mainLim = lims[int(maxlimInd[0])]
ax.set_xlim(mainLim)
ax.set_ylim(mainLim)
ax.set_zlim(mainLim)
plt.title('Initial Halo')
plt.savefig('init_halo')

plt.show()

ti = 0 #initial time
tf = 200 #final time, Myear
h  = 1 #timestep, Myear
l = 1e9 #box size, kpc
nSave = 1 #number of steps between each save
N = 3 #number of physical dimensions (plotting set for 3)
t1 = time.time()
bodlist, U_list, K_list = integrate(bodz,ti,tf,h,N,l,nSave)
t2 = time.time()
dt = t2-t1
print('integration finished after '+str(dt)+' seconds')

nP = len(bodlist[0])
pos1 = np.zeros((nP,N))
pos2 = np.zeros((nP,N))
vel2 = []


bodz1 = bodlist[0]
bodz2 = bodlist[-1]



for i in range(nP):
    for j in range(N):
        pos = bodz1[i].pos
        vel = np.sqrt(np.dot(bodz2[i].vel,bodz2[i].vel))
        pos1[i,j] = pos[j]
        vel2.append(vel)

print('!!!')
print(np.average(vel2))
print(min(vel2))
print(max(vel2))

pos1x = pos1[:,0]
pos1y = pos1[:,1]
pos1z = pos1[:,2]

pos1xR = pos1x[halo1_ind]
pos1yR = pos1y[halo1_ind]
pos1zR = pos1z[halo1_ind]

pos1xB = pos1x[halo2_ind]
pos1yB = pos1y[halo2_ind]
pos1zB = pos1z[halo2_ind]

# fig = plt.figure()
# ax = fig.add_subplot(111,projection='3d')
# ax.plot(pos1xR,pos1yR,pos1zR,'r.',markersize=2)
# ax.plot(pos1xB,pos1yB,pos1zB,'b.',markersize=2)
# plt.title('initial')
# ax.set_xlim(mainLim)
# ax.set_ylim(mainLim)
# ax.set_zlim(mainLim)

for i in range(nP):
    for j in range(N):
        pos = bodz2[i].pos
        pos2[i,j] = pos[j]

pos2x = pos2[:,0]
pos2y = pos2[:,1]
pos2z = pos2[:,2]

pos2xR = pos2x[halo1_ind]
pos2yR = pos2y[halo1_ind]
pos2zR = pos2z[halo1_ind]

pos2xB = pos2x[halo2_ind]
pos2yB = pos2y[halo2_ind]
pos2zB = pos2z[halo2_ind]

fig = plt.figure()
ax = fig.add_subplot(111,projection='3d')
ax.plot(pos2xR,pos2yR,pos2zR,'r.',markersize=2)
ax.plot(pos2xB,pos2yB,pos2zB,'b.',markersize=2)
ax.set_xlim(mainLim)
ax.set_ylim(mainLim)
ax.set_zlim(mainLim)
plt.title('Final Halo')
plt.savefig('final_halo')

Kvals = np.array(K_list)
Uvals = np.array(U_list)

Etot  = Kvals+Uvals
Etot0 = Etot[0]
Ediff = Etot - Etot0

virialE = 2*Kvals + Uvals

steps = np.arange(len(bodlist))
plt.figure()
plt.plot(steps,np.array(U_list)+np.array(K_list),'g',label='total energy')
plt.plot(steps,U_list,'b',label = 'potential energy')
plt.plot(steps,K_list,'r',label = 'kinetic energy')
plt.xlabel('timesteps')
plt.ylabel('energy (Mstar*(kpc/Gy)**2)')
plt.legend()
plt.title('All Energy')
plt.savefig('energy_all')

plt.figure()
plt.plot(steps,Ediff,'g',label='energy change')
plt.plot((0,steps[-1]),(0,0),'k:')
plt.xlabel('timesteps')
plt.ylabel('energy (Mstar*(kpc/Gy)**2)')
plt.legend()
plt.title('Total Energy Change')
plt.savefig('energy_total')

plt.figure()
plt.plot(steps,virialE,'b',label='2T+U')
plt.plot((0,steps[-1]),(0,0),'k:')
plt.xlabel('timesteps')
plt.ylabel('energy (Mstar*(kpc/Gy)**2)')
plt.legend()
plt.title('Virial Energy')
plt.savefig('energy_virial')
plt.show()

# format output in the way we need it for analysis
# graphs?? plt.plot(physics)
