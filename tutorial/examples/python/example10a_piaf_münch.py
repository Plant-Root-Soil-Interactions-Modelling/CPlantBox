import math
import os
import sys
sys.path.append("../../.."); sys.path.append("../../../src/python_modules")
sys.path.append("../../../tutorial/jupyter")
import datetime
import matplotlib.pylab as plt
import numpy as np
import timeit
import plotly
import math
import plotly.graph_objects as go
import xml.etree.ElementTree as ET
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import proj3d
from plotly.subplots import make_subplots

from vtk.util import numpy_support as VN
import vtk

import plantbox as pb
from phloem_flux import PhloemFluxPython
import vtk_plot as vp
#from rb_tools import * # CRootBox nodes Conversion tools
from CPlantBox_PiafMunch_new import * #adapted from xiaoran's CPlantBox_PiafMunch file

home_dir = os.getcwd()
dir_name = "/results"
dir_name2 = home_dir + dir_name
test = os.listdir(dir_name2)
for item in test:
    if item.endswith("10a.txt"):
        os.remove(os.path.join(dir_name2, item))
    if item.endswith("10a.vtk"):
        os.remove(os.path.join(dir_name2, item))
    if item.endswith("10a.vtp"):
        os.remove(os.path.join(dir_name2, item))


name = "oneroot_10a" # parameter name
time = 1 
plant = pb.MappedPlant()
path = "../../../modelparameter/plant/"
plant.readParameters(path + name + ".xml")
plant.initialize()
plant.simulate(time, False)
phl = PhloemFluxPython(plant)
ana = pb.SegmentAnalyser(phl.rs)
ana.write("results/example10a.vtp" )

write_PiafMunch_parameter(phl) 
 




