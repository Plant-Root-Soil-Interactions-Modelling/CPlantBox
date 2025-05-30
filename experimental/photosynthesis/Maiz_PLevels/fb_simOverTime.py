import sys;

CPBdir = "../../.."
sys.path.append(CPBdir + "/src");
sys.path.append(CPBdir);
sys.path.append("../../..");
sys.path.append("..");
sys.path.append(CPBdir + "/src/python_modules");
sys.path.append("../build-cmake/cpp/python_binding/")  # dumux python binding
sys.path.append("../../build-cmake/cpp/python_binding/")
sys.path.append("../modules/")  # python wrappers
sys.path.append("../../experimental/photosynthesis/")

import visualisation.vtk_plot as vp
import importlib
import pandas as pd
import plantbox as pb
from functional.phloem_flux import PhloemFluxPython  # Python hybrid solver
import math
import os
import numpy as np
from datetime import datetime, timedelta

from help_MaizPlevels import setKrKx_xylem, sinusoidal, theta2H, setKrKx_phloem
# reload(helpuqrMasterCopy1)
import numpy as np




def initPlant():
    simMax = 28
    # plant system
    pl = pb.MappedPlant(seednum=2)  # pb.MappedRootSystem() #pb.MappedPlant()
    # pl2 = pb.MappedPlant(seednum = 2) #pb.MappedRootSystem() #pb.MappedPlant()
    path = CPBdir + "/modelparameter/structural/plant/"
    for Plevel in [0,3]:
        name = "P"+str(Plevel)+"_plant"  # "Triticum_aestivum_adapted_2021"#
        pl.readParameters( name + ".xml")

        if name == 'P0_plant':
            for p in pl.getOrganRandomParameter(pb.leaf):
                p.la,  p.lmax = 38.41053981, 38.41053981
                p.areaMax = 54.45388021  # cm2, area reached when length = lmax
                NLeaf = 100  # N is rather high for testing
                phi = np.array([-90,-80, -45, 0., 45, 90]) / 180. * np.pi    
                l = np.array([38.41053981,1 ,1, 0.3, 1, 38.41053981]) #distance from leaf center
                p.tropismT = 1 # 6: Anti-gravitropism to gravitropism
                #p.tropismN = 5
                p.tropismS = 0.05
                p.tropismAge = 5 #< age at which tropism switch occures, only used if p.tropismT = 6
                p.createLeafRadialGeometry(phi,l,NLeaf)
            for p in pl.getOrganRandomParameter(pb.stem):
                p.r = 0.758517633
                p.lmax = (simMax-7)*0.758517633 
                # p.lmax=200
        if name == 'P1_plant':
            for p in pl.getOrganRandomParameter(pb.leaf):
                p.lb =  0 # length of leaf stem
                p.la,  p.lmax = 42.60617256, 42.60617256
                p.areaMax = 66.69532685  # cm2, area reached when length = lmax
                NLeaf = 100  # N is rather high for testing
                phi = np.array([-90,-80, -45, 0., 45, 90]) / 180. * np.pi
                l = np.array([42.60617256,1 ,1, 0.3, 1, 42.60617256]) #distance from leaf center
                p.tropismT = 1 # 6: Anti-gravitropism to gravitropism
                #p.tropismN = 5
                p.tropismS = 0.05
                p.tropismAge = 5 #< age at which tropism switch occures, only used if p.tropismT = 6
                p.createLeafRadialGeometry(phi, l, NLeaf)
            for p in pl.getOrganRandomParameter(pb.stem):
                r= 0.91546738
                p.r = r
                # p.lmax = (simMax-7)*0.91546738  
        if name == 'P2_plant':
            for p in pl.getOrganRandomParameter(pb.leaf):
                p.lb =  0 # length of leaf stem
                p.la,  p.lmax = 52.23664394, 52.23664394
                p.areaMax = 80.68274258  # cm2, area reached when length = lmax
                NLeaf = 100  # N is rather high for testing
                phi = np.array([-90,-80, -45, 0., 45, 90]) / 180. * np.pi
                l = np.array([52.23664394,1 ,1, 0.3, 1, 52.23664394]) #distance from leaf center
                p.tropismT = 1 # 6: Anti-gravitropism to gravitropism
                #p.tropismN = 5
                p.tropismS = 0.05
                p.tropismAge = 5 #< age at which tropism switch occures, only used if p.tropismT = 6
                p.createLeafRadialGeometry(phi, l, NLeaf)
            for p in pl.getOrganRandomParameter(pb.stem):
                r= 1.000613891
                p.r = r
                p.lmax = (simMax-7)*1.000613891    

        if name == 'P3_plant':
            for p in pl.getOrganRandomParameter(pb.leaf):
                p.lb =  0 # length of leaf stem
                p.la,  p.lmax = 49.12433414, 49.12433414
                p.areaMax = 71.95670914  # cm2, area reached when length = lmax
                NLeaf = 100  # N is rather high for testing
                phi = np.array([-90,-80, -45, 0., 45, 90]) / 180. * np.pi
                l = np.array([49.12433414,1 ,1, 0.3, 1, 49.12433414]) #distance from leaf center
                p.tropismT = 1 # 6: Anti-gravitropism to gravitropism
                p.tropismN = 5
                p.tropismS = 0.05
                p.tropismAge = 5 #< age at which tropism switch occures, only used if p.tropismT = 6
                p.createLeafRadialGeometry(phi, l, NLeaf)

            for p in pl.getOrganRandomParameter(pb.stem):
                r= 1.128705967
                p.r = r
                # p.lmax = (simMax-7)*1.128705967   
        pl.writeParameters(name+ "new.xml")
initPlant()