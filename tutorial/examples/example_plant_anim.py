"""plant example"""
import sys ;sys.path.append("../../.."); sys.path.append("../../../src/python_modules")
import numpy as np;

import plantbox as pb
import vtk_plot as vp
from astropy.units import dd

plant = pb.Plant()

# Open plant and root parameter from a file
path = "../../../modelparameter/plant/"
name = "0"  # CPlantBox_test_leaf_tree00

# LEAFS smallPlant_mgiraud "manyleaves"
# NO LEAFS "CPlantBox_test_leaf_tree22"  # "chicon_entire"  # "Anagallis_femina_leaf_shape"  # "Anagallis_femina_Leitner_2010"

# BREAKS MY COMPUTER Swiss_chard

plant.readParameters(path + name + ".xml")

# print radii
print("\nleafs")
for p in plant.getOrganRandomParameter(pb.leaf):
    p.a = 0.05
    p.a_s = 0
    if (p.subType > 0): 
        print(p.subType, "radius", p.a, "lmax", p.lmax, p.ln, p.lb, p.successor, p.nob())    
        
        if (p.subType == 2): 
            p.tropismT = 1
            p.tropismN = 5
            p.tropismS = 0.01  # 0.3
            
        if (p.subType == 3): 
            # print(p)
            
            p.la, p.lb, p.lmax, p.ln, = 3.5, 1., 7.5, 3  
            p.areaMax = 10  # cm2
            phi = np.array([-90, -45, 0., 45, 90]) / 180. * np.pi
            l = np.array([3, 2.2, 1.7, 2, 3.5])
            N = 101  # N is rather high for testing
            p.createLeafRadialGeometry(phi, l, N)
      
#             lrp = p
#             p.areaMax = 20  # cm2            
#             lrp.la, lrp.lb, lrp.lmax, lrp.ln, lrp.r, lrp.dx = 5, 1, 11, 5, 1, 0.1
#             phi = np.array([-90., -67.5, -45, -22.5, 0, 22.5, 45, 67.5, 90]) / 180. * np.pi
#             l = np.array([5., 1, 5, 1, 5, 1, 5, 1, 5])
#             assert(l.shape == phi.shape)
#             N = 500  # N is rather high for testing
#             lrp.createLeafRadialGeometry(phi, l, N)               
#             p = lrp
            
            p.tropismT = 1
            p.tropismN = 5
            p.tropismS = 0.1  # 0.3
   
        else:
            p.a = p.a * 3
    
print("\nstem")
for p in plant.getOrganRandomParameter(pb.stem):
    print(p.subType, "radius", p.a, "lmax", p.lmax, p.ln, p.lb, p.successor, p.nob())  
    if (p.subType == 1): 
        p.a = p.a / 2
        p.a = 0.2
        p.a_s = 0     
        p.lmax = 30
        p.la = 3
        p.ln = 1

print("roots")
for p in plant.getOrganRandomParameter(pb.root):
    if (p.subType > 0): 
        print(p.subType, p.a, p.successor)
        if p.subType == 1:  # taproot
            p.theta = 0.
        p.a = 0.05
        p.a_s = 0

# soil = pb.SDF_PlantContainer(1.e6, 1.e6, 1.e6, False)
# plant.setGeometry(soil)

# increase resolution
for p in plant.getOrganRandomParameter(pb.root):
    p.dx = 0.2

# Initialize
plant.initialize()

dt = 0.1
N = 400
min_ = np.array([-20, -20, -50])
max_ = np.array([20, 20, 30.])

# test = plant.getOrgans(pb.leaf)
# print("test")

anim = vp.AnimateRoots(plant)
anim.min = min_
anim.max = max_
anim.res = [1, 1, 1]
anim.file = "results/example_plant"
anim.avi_name = "results/example_"
anim.plant = True
anim.start()

for i in range(0, N):

    plant.simulate(dt, False)
    anim.root_name = "organType"
    anim.update()

