"""scales insertion angle"""
import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import visualisation.vtk_plot as vp
import numpy as np

plant = pb.Plant()
path = "../../modelparameter/structural/rootsystem/"
name = "Glycine_max_Moraes2020"
plant.readParameters(path + name + ".xml")

box = pb.SDF_PlantBox(10, 10, 30)  # nutrient rich patch  # |\label{l34:patch2}|
patch = pb.SDF_RotateTranslate(box, pb.Vector3d(-5, 0., -10))

max_ = 1.  # maximal |\label{l34:rate2_start}|
min_ = 0.1  # minimal
slope = 1.  # [cm] linear gradient between min and max
soilprop = pb.SoilLookUpSDF(patch, max_, min_, slope)  # |\label{l34:rate2}|

for organ_type in [pb.root, pb.stem, pb.leaf]:
    for p in plant.getOrganRandomParameter(organ_type):
        if p.subType > 2:
            p.dx = 0.25  # adjust resolution
            p.f_sa = soilprop  # Scale insertion angle
            p.lmax = 2 * p.lmax  # make second order laterals longer

plant.initialize()
simtime = 15.
dt = 1.
for i in range(0, round(simtime / dt)):
    plant.simulate(dt)

ana = pb.SegmentAnalyser(plant)  # |\label{l34:analysis2}|
ana.crop(patch)
lm_theta = np.mean(ana.getParameter("theta"))
ana = pb.SegmentAnalyser(plant)
ana.crop(pb.SDF_Complement(patch))
rm_theta = np.mean(ana.getParameter("theta"))
print('\nMean insertion angle within patch {:g} degrees'.format(lm_theta))
print('Mean insertion angle outside patch {:g} degrees\n'.format(rm_theta))

plant.write("results/example_3_4b.vtp")
vp.plot_roots_and_container(plant, patch, "theta")
