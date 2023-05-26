import sys; sys.path.append("../.."); sys.path.append("../../src/")
import numpy as np
import matplotlib.pyplot as plt

import plantbox as pb
import visualisation.vtk_plot as vp


pl = pb.MappedPlant()  
CPBdir = "../.."
path = CPBdir + "/modelparameter/structural/plant/"
name = "pseudoStem"  
pl.readParameters( path+name + ".xml")
pl.initialize(verbose=False) 
N = 100
#for i in range(N):
#    pl.simulate(1)
#    #ana = pb.SegmentAnalyser(pl.mappedSegments())
#    #ana.write("pseudoStem_"+str(i)+".vtp", ["organType","subType"]) 
    
    
pl.simulate(100)
vp.plot_plant(pl, "organType")