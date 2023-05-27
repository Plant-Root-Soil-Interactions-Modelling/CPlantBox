import sys; sys.path.append("../.."); sys.path.append("../../src/")
import numpy as np
import matplotlib.pyplot as plt

import plantbox as pb
import visualisation.vtk_plot as vp


pl = pb.MappedPlant()  
mpl = pb.MappedPlant()  
CPBdir = "../.."
path = CPBdir + "/modelparameter/structural/plant/"
name = "Triticum_aestivum_test_2021"  
pl.readParameters( path+name + ".xml")
pl.initialize(verbose=False) 
mpl.readParameters( path+name + ".xml")
mpl.initialize(verbose=False) 
N = 50
mpl.simulate(5, False)
for i in range(N):
    #pl.simulate(1)
    mpl.simulate(1, False)
    organTypes = np.array(mpl.organTypes)
    subTypes   = np.array(mpl.subTypes)
    
    #mana = pb.SegmentAnalyser(mpl.mappedSegments())
    #mana.write("pseudoStem_m"+str(i)+".vtp", ["organType","subType"]) 
    leafes = mpl.getOrgans(pb.leaf)
    print(organTypes, subTypes)
    if len(leafes) > 0:
        vp.plot_plant(mpl,p_name = [ "organType","subType"],
                            vals =[ organTypes,subTypes], 
                            filename = "pseudoStem_m"+str(i),
                            range_ = [0,5000])
#ana = pb.SegmentAnalyser(pl.mappedSegments())
#ana = pb.SegmentAnalyser(pl)

#ana.addData("ot", x)
#ana.write("pseudoStem_"+str(i)+".vtp", ["organType","subType"]) 

pl.simulate(N, False)   
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