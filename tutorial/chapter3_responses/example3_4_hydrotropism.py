"""hydrotropism"""

import plantbox as pb
import plantbox.visualisation.vtk_plot as vp
import numpy as np

plant = pb.Plant()
path = "../../modelparameter/structural/rootsystem/"
filename = "Glycine_max_Moraes2020"
plant.readParameters(path + filename + ".xml")

box = pb.SDF_PlantBox(4, 4, 4)  # nutrient rich patch  # |\label{l36:patch2}|
patch = pb.SDF_RotateTranslate(box, pb.Vector3d(-20, -20, -30))

max_ = -300  # maximal
min_ = -10000  # minimal
slope = 0.01  # linear gradient between min and max (cm)
soilprop = pb.SoilLookUpSDF(patch, max_, min_, slope) 

for p in plant.getOrganRandomParameter(pb.root):
    p.tropismT = pb.TropismType.hydro
    p.tropismN = 2                          # number of trials per segment
    p.tropismS = 0.1           # sensitivity per root type  # |\label{l36:sigma}|

plant.setSoil(soilprop)                # |\label{l36:setsoil}|

plant.initialize()  # |\label{l36:loop}|
sim_time = 15.0
dt = 1.0
for i in range(0, round(sim_time / dt)):
    plant.simulate(dt)  # |\label{l36:loop_end}|


vp.plot_roots_and_container(plant, patch, "subType")
