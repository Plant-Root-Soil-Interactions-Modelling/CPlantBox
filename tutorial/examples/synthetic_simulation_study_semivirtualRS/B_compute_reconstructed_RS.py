"""small example"""
import sys
sys.path.append("../../../")
sys.path.append("../../../src/");
sys.path.append("../../../gui/viewer/");
import plantbox as pb
from functional.xylem_flux import XylemFluxPython 
import visualisation.vtk_plot as vp
import visualisation.vtk_tools as vt
from viewer_data import ViewerDataModel
from visualisation.vtk_tools import *
import rsml.rsml_writer
import rsml.rsml_reader as rsml
import numpy as np
import matplotlib.pyplot as plt
import os
from random import randint
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

        #query nodes, segs 
        ana = pb.SegmentAnalyser(rs)
        #filter out all 3rd order roots 
        ana.filter("subType", 1, 2)
        nodes = np_convert(ana.nodes)
        nodes = np.around(nodes, decimals=2)
        segs = np_convert(ana.segments)

        #store all original measures 
        radius = np.array(ana.getParameter("radius"))
        length = np.array(ana.getParameter("length"))
        cts = np.array(ana.getParameter("creationTime"))
        types = np.array(ana.getParameter("subType"), dtype = int)
        
        radius = radius.tolist()
        nodes_ = nodes.tolist()
        cts = cts.tolist()
        segs = segs.tolist()
        subtype = types.tolist()
        length = length.tolist()

        nodes = [pb.Vector3d(n[0], n[1], n[2]) for n in nodes]
        segs = [pb.Vector2i(s[0], s[1]) for s in segs]
        
        ana = pb.SegmentAnalyser(nodes, segs, cts, radius)
        ana.addData("radius", radius)
        ana.addData("subType", subtype)
        ana.addData("creationTime", cts)
        ana.addData("length", length)
        fname = "results/REC_"+namep[0]+'_day'+str(ii)
        ana.write(fname+".vtp", ["radius", "subType", "creationTime", "length"] )

        #write rsml
        pd = vt.read_vtp(fname +'.vtp')
        vt.write_rsml(fname + '.rsml', pd,0,[],[1])
        print('done')
        



 
