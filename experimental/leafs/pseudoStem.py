import sys; sys.path.append("../.."); sys.path.append("../../src/")
import numpy as np
import matplotlib.pyplot as plt

import plantbox as pb
import visualisation.vtk_plot as vp


pl = pb.MappedPlant()  
mpl = pb.MappedPlant()  
CPBdir = "../.."
path = CPBdir + "/modelparameter/structural/plant/"
name = "pseudoStem"  
pl.readParameters( path+name + ".xml")
pl.initialize(verbose=False) 
mpl.readParameters( path+name + ".xml")
mpl.initialize(verbose=False) 
N = 100
for i in range(N):
    #pl.simulate(1)
    mpl.simulate(1, True)
    mana = pb.SegmentAnalyser(mpl.mappedSegments())
    mana.write("pseudoStem_m"+str(i)+".vtp", ["organType","subType"]) 
#ana = pb.SegmentAnalyser(pl.mappedSegments())
#ana = pb.SegmentAnalyser(pl)

#ana.addData("ot", x)
#ana.write("pseudoStem_"+str(i)+".vtp", ["organType","subType"]) 

pl.simulate(N, True)   
ana = pb.SegmentAnalyser(pl.mappedSegments())
ana.write("pseudoStem_"+str(i)+".vtp", ["organType","subType"]) 
leaves = mpl.getOrgans(pb.leaf)   
if False:
    for leaf in leaves:
        print(np.array([np.array(leaf.getNode(i)) for i in range(leaf.getNumberOfNodes())]))
    leaves = mpl.getOrgans(pb.stem)   
    for leaf in leaves:
        print(np.array([np.array(leaf.getNode(i)) for i in range(leaf.getNumberOfNodes())]))
#print(ana.data)
#print(np.array([np.array(i) for i in pl.segments]))
#print(np.array([np.array(i) for i in mpl.segments]))
#pl.simulate(100)
#vp.plot_plant(pl, "organType")
#vp.plot_plant(mpl, "organType")