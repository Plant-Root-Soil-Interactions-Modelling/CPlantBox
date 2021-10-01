

import sys; sys.path.append(".."); sys.path.append("../src/python_modules");sys.path.append("../src")
import plantbox as pb
import vtk_plot as vp

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
    if item.endswith("NodalGrowth.vtp"):
        os.remove(os.path.join(item))

    
    
######################
#
# plant
#
####################### 
pl = pb.MappedPlant(seednum = 1)#set seed
path = "../modelparameter/plant/" #"../../../modelparameter/rootsystem/" 
name = "4testrel"

pl.readParameters(path + name + ".xml")
dt =1
steps =35
pl.initialize(verbose = True)

for step in range(steps):
    print("\n\n\nstep n°", step)
    pl.simulate(dt, True)
    ana = pb.SegmentAnalyser(pl)
    ana.write("%s_NodalGrowth.vtp" %(step))
