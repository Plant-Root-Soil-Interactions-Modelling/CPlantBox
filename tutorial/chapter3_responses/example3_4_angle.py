"""scales insertion angle"""

import numpy as np

import plantbox as pb
import plantbox.visualisation.vtk_plot as vp

plant = pb.Plant()
path = "../../modelparameter/structural/rootsystem/"
name = "Glycine_max_Moraes2020"
plant.readParameters(path + name + ".xml")

box = pb.SDF_PlantBox(10, 10, 30)  # nutrient rich patch  # |\label{l34:patch2}|
patch = pb.SDF_RotateTranslate(box, pb.Vector3d(-5, 0.0, -10))

max_ = 1.0  # maximal
min_ = 0.1  # minimal
slope = 1.0  # linear gradient between min and max (cm) 
soilprop = pb.SoilLookUpSDF(patch, max_, min_, slope)  # |\label{l34:rate2}|

for organ_type in [pb.root, pb.stem, pb.leaf]:  # |\label{l34:for_start}|
    for p in plant.getOrganRandomParameter(organ_type):
        if p.subType > 2:
            p.dx = 0.25  # adjust resolution (cm)
            p.f_sa = soilprop  # scale insertion angle
            p.lmax = 2 * p.lmax  # increase higher order length

plant.initialize()
simtime = 15.0
dt = 1.0
for i in range(0, round(simtime / dt)):
    plant.simulate(dt)

ana = pb.SegmentAnalyser(plant)  # |\label{l34:analysis2}|
ana.crop(patch)
lm_theta = np.mean(ana.getParameter("theta_deg"))
ana = pb.SegmentAnalyser(plant)
ana.crop(pb.SDF_Complement(patch))
rm_theta = np.mean(ana.getParameter("theta_deg"))
print(f"\nMean insertion angle within patch {lm_theta:g} degrees")
print(f"Mean insertion angle outside patch {rm_theta:g} degrees\n")

plant.write("results/example_3_4b.vtp")
vp.plot_roots_and_container(plant, patch, "theta_deg")
