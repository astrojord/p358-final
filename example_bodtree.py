import numpy as np
from BHA import Body
from BHA import Node
from integrator import grav_accelerate

tree = Node(pos = [0,0,0], length = 10)
for i in range(100000):
    tree.fit([10*np.random.random(3) - 5., np.random.random()])
tree.calculate_coms()
bod = Body((2,3,4),6)
mrrogers = tree.neighbors(bod)
print(grav_accelerate(bod,tree))

#mr rogers is a list of all neighbors
#where the neighbor data are tuples with a position array and mass float
