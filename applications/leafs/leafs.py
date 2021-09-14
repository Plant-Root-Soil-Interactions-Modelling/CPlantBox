import sys; sys.path.append(".."); sys.path.append("../../src/python_modules"); sys.path.append("../../")
import numpy as np
import matplotlib.pyplot as plt

import plantbox as pb
import vtk_plot as vp

plant = pb.Plant()  

seed_rp = pb.SeedRandomParameter(plant)
seed_rp.subType = 0
plant.setOrganRandomParameter(seed_rp)

root_rp = pb.RootRandomParameter(plant)
root_rp.subType = 1
root_rp.lmax = 1
root_rp.r = 0.1
root_rp.theta = 0
plant.setOrganRandomParameter(root_rp)

stem_rp = pb.StemRandomParameter(plant)
stem_rp.subType = 1
stem_rp.la = 1.
stem_rp.lb = 5
stem_rp.lmax = 7.5
stem_rp.ln = 1
stem_rp.theta = 0
stem_rp.successor = [2]
stem_rp.successorP = [1]

plant.setOrganRandomParameter(stem_rp)
stem_rp = pb.StemRandomParameter(plant)
stem_rp.subType = 2  # dummy 
stem_rp.la = 0
stem_rp.lb = 0
stem_rp.lmax = 0
stem_rp.ln = 0
plant.setOrganRandomParameter(stem_rp)

leaf_rp = pb.LeafRandomParameter(plant)
leaf_rp.subType = 2
leaf_rp.lmax = 10

leaf_rp.tropismS = 0.
     
leaf_rp.theta = 20. / np.pi

# Mint
# leaf_rp.la, leaf_rp.lb, leaf_rp.lmax, leaf_rp.ln, leaf_rp.r, leaf_rp.dx = 3.5, 1., 7.5, 3, 1, 0.5  
# phi = np.array([-90, -45, 0., 45, 90]) / 180. * np.pi
# l = np.array([3, 2.2, 1.7, 2, 3.5])
# N = 30  # N is rather high for testing
# leaf_rp.createLeafRadialGeometry(phi, l, N)    

# Pasley (non convex - TODO)
# leaf_rp.la, leaf_rp.lb, leaf_rp.lmax, leaf_rp.ln, leaf_rp.r, leaf_rp.dx = 5.4, 4.3, 5.4 + 4.3 + 0.8, 0.8, 1, 0.5
# phi = np.array([-90, -45, -25, -15, 0, 10, 20, 30, 40, 45, 55, 60, 65, 90]) / 180. * np.pi
# l = np.array([0.8, 0.5, 1.5, 2, 2.6, 2, 2.6, 2.7, 3.1, 3.3, 2.6, 3.3, 3.8, 5.4])
# assert(l.shape == phi.shape)
# N = 200  # N is rather high for testing
# leaf_rp.createLeafRadialGeometry(phi, l, N)    

leaf_rp.la, leaf_rp.lb, leaf_rp.lmax, leaf_rp.ln, leaf_rp.r, leaf_rp.dx = 5, 1, 11, 5, 1, 0.25
phi = np.array([-90., -67.5, -45, -22.5, 0, 22.5, 45, 67.5, 90]) / 180. * np.pi
l = np.array([5., 2, 5, 2, 5, 2, 5, 2, 5])
assert(l.shape == phi.shape)
N = 200  # N is rather high for testing
leaf_rp.createLeafRadialGeometry(phi, l, N)   

leaf_rp.areaMax = 50  # cm2

plant.setOrganRandomParameter(leaf_rp)        

plant.initialize()
plant.simulate(100)

for o in plant.getOrgans():
    print(o)

vp.plot_plant(plant, "creationTime")

print("fin")

