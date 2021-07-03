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

""" Parameters """
kz = 4.32e-1  # axial conductivity [cm^3/day] 
kr = 1.728e-4  # radial conductivity of roots [1/day]
kr_stem = 1.e-20  # radial conductivity of stem  [1/day], set to almost 0
gmax = 0.004 #  cm3/day radial conductivity between xylem and guard cell
p_s = -200  # static water potential (saturation) 33kPa in cm
#p_g = -2000 # water potential of the guard cell
RH = 0.5 # relative humidity
TairC = 20
p_a =  -1000  #default outer water potential 
simtime = 14.0  # [day] for task b
k_soil = []
Q = 900e-6 # mol quanta m-2 s-1 light, example from leuning1995
cs = 350e-6 #co2 paartial pressure at leaf surface (mol mol-1)
TairK = TairC + 273.15
p_linit = -27000
ci_init = cs * 0.3

es = 0.61078 * math.exp(17.27 * TairC / (TairC + 237.3)) 
ea = es * RH 
VPD = es - ea 


name = "Triticum_aestivum_adapted_2021"#"morning_glory_7m"#"Triticum_aestivum_adapted_2021"#"morning_glory_7m"#"oneroot_10a" # parameter name
time = 35
plant = pb.MappedPlant()
path = "../../../modelparameter/plant/"
plant.readParameters(path + name + ".xml")
plant.initialize()
plant.simulate(time, False)
phl = PhloemFluxPython(plant)
ana = pb.SegmentAnalyser(phl.rs)
ana.write("results/example10a.vtp" )

nodes = phl.get_nodes()
tiproots, tipstem, tipleaf = phl.get_organ_nodes_tips() #end node of end segment of each organ
node_tips = np.concatenate((tiproots, tipstem, tipleaf))
tiproots, tipstem, tipleaf = phl.get_organ_segments_tips() #end segment of each organ
seg_tips = np.concatenate((tiproots, tipstem, tipleaf))
initphi=[]
newinit = False
if newinit:
    phl.setKr([[kr],[kr_stem],[gmax]]) #att: check units
    phl.setKx([[kz]])
    phl.airPressure = p_a

    phl.seg_ind = seg_tips # segment indices for Neumann b.c.
    phl.node_ind = node_tips

    rx = phl.solve_leuning(sim_time = simtime,sxx=[p_s], cells = True, Qlight = Q,VPD = VPD,
                    Tl = TairK,p_linit = p_linit,ci_init = ci_init,cs=cs, soil_k = [], log = True)
                    
    minPhi = 2e-4 
    maxPhi = 8e-4

    #link init C val to rx => makes it linked to the distance to the leaf
    coef = (minPhi-maxPhi)/(max(rx)-min(rx))
    param = maxPhi - min(rx)*coef

    initphi = rx*coef + param
    #print(initphi)
    ana = pb.SegmentAnalyser(phl.rs)
    ana.addData("initPhi", initphi) 
    ana.write("results/example10a.vtp" , [ "phi"])
nopsi=True
name="Triticum_a_"
write_PiafMunch_parameter(phl, initphi, name = name, nopsi=nopsi) 
 


