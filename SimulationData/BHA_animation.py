import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from BHA import Body
from BHA import Node
import matplotlib as mpl
from matplotlib import animation
import ffmpeg

mpl.rcParams['lines.linewidth'] = 0.1

data = np.loadtxt('snapshot100.txt')

posfull = data[:,0:3]
pos = []
for i in range(data.shape[0]):
    if i%50 == 0:
        #pos.append(posfull[i,:])
        p = (np.random.rand(1,3)-0.5)*(i+1)
        print(p)
        pos.append(p[0])

pos = []

for i in range(100):
    a = (np.random.rand(1,3)-0.5)*(10*(i+1)**0.1)
    b = (np.random.rand(1,3))*(7*(i+1)**0.1)-[12,34,-19]
    pos.append(a[0])
    if i%3 == 0:
        pos.append(b[0])


pos = np.array(pos)
print(pos.shape)
fig = plt.figure()
ax = fig.add_subplot(111,projection='3d')
ax.plot(pos[:,0],pos[:,1],pos[:,2],'r.',markersize=2)

objlist = []
for i in range(pos.shape[0]):
    obj = Body(pos[i,:],1,[0,0,0],[0,0,0],1)
    objlist.append(obj)

tree = Node(pos = [0, 0, 0], length = 75)
tree.fit(objlist)
tree.calculate_coms()

L = tree.length
mainLim = (-L/2,L/2)

def gridplot(node):
    children = node.children
    for child in children:
        #print('child!')
        center = child.pos
        cx,cy,cz = center[0],center[1],center[2]
        l = child.length

        # fig = plt.figure()
        # ax = fig.add_subplot(111,projection='3d')
        # ax.plot(pos[:,0],pos[:,1],pos[:,2],'r.',markersize=2)

        ax.plot((cx+l/2,cx+l/2),(cy+l/2,cy+l/2),(cz-l/2,cz+l/2),'k')
        ax.plot((cx+l/2,cx+l/2),(cy-l/2,cy-l/2),(cz-l/2,cz+l/2),'k')
        ax.plot((cx-l/2,cx-l/2),(cy+l/2,cy+l/2),(cz-l/2,cz+l/2),'k')
        ax.plot((cx-l/2,cx-l/2),(cy-l/2,cy-l/2),(cz-l/2,cz+l/2),'k')

        ax.plot((cx+l/2,cx+l/2),(cy-l/2,cy+l/2),(cz+l/2,cz+l/2),'k')
        ax.plot((cx+l/2,cx+l/2),(cy-l/2,cy+l/2),(cz-l/2,cz-l/2),'k')
        ax.plot((cx-l/2,cx-l/2),(cy-l/2,cy+l/2),(cz+l/2,cz+l/2),'k')
        ax.plot((cx-l/2,cx-l/2),(cy-l/2,cy+l/2),(cz-l/2,cz-l/2),'k')

        ax.plot((cx-l/2,cx+l/2),(cy+l/2,cy+l/2),(cz+l/2,cz+l/2),'k')
        ax.plot((cx-l/2,cx+l/2),(cy+l/2,cy+l/2),(cz-l/2,cz-l/2),'k')
        ax.plot((cx-l/2,cx+l/2),(cy-l/2,cy-l/2),(cz+l/2,cz+l/2),'k')
        ax.plot((cx-l/2,cx+l/2),(cy-l/2,cy-l/2),(cz-l/2,cz-l/2),'k')

        ax.set_xlim(mainLim)
        ax.set_ylim(mainLim)
        ax.set_zlim(mainLim)

        #plt.show()
        #print(child.type)
        if child.type == 'INTERNAL' and l > 0.1*L:
            #print('rec')
            #print(child.type)
            gridplot(child)

#gridplot(tree)
for l in fig.gca().lines:
    l.set_alpha(1)

#Writer = animation.writers['ffmpeg']
#writer = Writer(fps=30, metadata=dict(artist='Me'), bitrate=1800)

def init():
    #ax.scatter(xx, yy, zz, marker='o', s=20, c="goldenrod", alpha=0.6)
    #gridplot(tree)
    ax
    return fig,

def animate(i):
    ax.view_init(elev=10., azim=i)
    return fig,

gridplot(tree)
# Animate
anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=360, interval=20, blit=True)
# Save
anim.save('BHA_ex2.html', fps=30, extra_args=['-vcodec', 'libx264'],)#writer=writer)

plt.show()
