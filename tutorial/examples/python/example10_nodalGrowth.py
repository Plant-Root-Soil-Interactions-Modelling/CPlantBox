

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

home_dir = os.getcwd()
dir_name = "/results"
dir_name2 = home_dir + dir_name
test = os.listdir(dir_name2)
for item in test:
    if item.endswith("10.txt"):
        os.remove(os.path.join(dir_name2, item))
    if item.endswith("10.vtk"):
        os.remove(os.path.join(dir_name2, item))
    if item.endswith("10.vtp"):
        os.remove(os.path.join(dir_name2, item))

    
    
######################
#
# plant
#
####################### 
pl = pb.MappedPlant(seednum = 1)#set seed
path = "../../../modelparameter/plant/" #"../../../modelparameter/rootsystem/" 
name = "4testrel"

pl.readParameters(path + name + ".xml")
dt =1
steps =35
pl.initialize(verbose = True)

for step in range(steps):
    print("\n\n\nstep n°", step)
    pl.simulate(dt, True)
    ana = pb.SegmentAnalyser(pl)
    ana.write("results/%s_example_10.vtp" %(step))

