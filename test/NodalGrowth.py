

import sys; sys.path.append(".."); sys.path.append("../src/python_modules");sys.path.append("../src")
#sys.path.append("/home/rbtlm2004/CPlantBoxinline/src");sys.path.append("/home/rbtlm2004/CPlantBoxinline/src/python_modules")
#sys.path.append("/home/rbtlm2004/miniconda3/envs/MYFIPYENV/lib/python3.8/site-packages/fipy_changed")

import plantbox as pb
import vtk_plot as vp
#from phloem_flux import PhloemFluxPython, GetTime
import math
import os
import pstats
from io import StringIO
from datetime import datetime, timedelta
import numpy as np
import random
       
class NullIO(StringIO):
    def write(self, txt):
       pass

test = os.listdir(os.getcwd())       
for item in test:
    if item.endswith("10b.txt"):
        os.remove(os.path.join(item))
    if item.endswith("10b.vtk"):
        os.remove(os.path.join( item))
    if item.endswith("10b.vtp"):
        os.remove(os.path.join(item))

    
    
######################
#
# plant
#
####################### 
#random.seed(10)
pl = pb.MappedPlant(seednum = 1)
#take care of MappedPlant afterwards
path = "../modelparameter/plant/" #"../../../modelparameter/rootsystem/" 
name = "4testrel"# "Triticum_aestivum_adapted_2021"# "Triticum_aestivum_adapted_2021"#"22ss"#"small_2020"
pl.readParameters(path + name + ".xml")
dt =1
steps =35
pl.initialize(verbose = True)

for step in range(steps):
    #print("\n\n\nstep n°", step)
    pl.simulate(dt, False)      
     
    ana = pb.SegmentAnalyser(pl)
    ana.write("%s_example10b.vtp" %(step))

