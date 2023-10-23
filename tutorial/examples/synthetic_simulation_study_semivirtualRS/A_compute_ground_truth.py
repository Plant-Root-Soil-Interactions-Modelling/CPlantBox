"""small example"""
import sys
sys.path.append("../../../")
sys.path.append("../../../src/");
sys.path.append("../../../gui/viewer/");
import plantbox as pb
from functional.xylem_flux import XylemFluxPython 
import visualisation.vtk_plot as vp
from viewer_data import ViewerDataModel
from visualisation.vtk_tools import *
import rsml.rsml_writer
import rsml.rsml_reader as rsml
import numpy as np
import matplotlib.pyplot as plt
import os
#########################################################
rs = pb.MappedRootSystem()
path = "../../../modelparameter/structural/rootsystem/"

name = "Faba_synMRI"
namep = name.split('_')
rs.readParameters(path + name + ".xml")
width = 2.8
depth = 19
RSage = 10 

#general info
bigcyl = pb.SDF_PlantContainer(width-0.05, width-0.05, depth-0.05, False) #depth -0.05

#initialize and simulate
rs.setGeometry(bigcyl)
rs.setSeed(1) #use the same seed for all simulations
rs.initialize()

for ii in range(0,RSage+1): 
    simtime = 1
    rs.simulate(simtime, False)

    
    if int(ii) == RSage:
        fname = "results/GT_"+namep[0]+'_day'+str(ii)
        rs.write(fname+".vtp")
        rs.write(fname+".rsml")



 
