"""increase axial resolution (e.g. for animation)"""
import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import visualisation.vtk_plot as vp

plant = pb.Plant()
path = "../../modelparameter/structural/plant/"
name = "fspm2023"  # "hello_world"
plant.readParameters(path + name + ".xml")

# Parameters for animation
sim_time = 30  # days
fps = 30  # frames per second
anim_time = 5  # seconds
N = fps * anim_time  # [1]
dt = sim_time / N  # days
filename = "animate"

# Modify axial resolution of all relevant organs
for organ_type in [pb.root, pb.stem, pb.leaf]:
    for p in plant.getOrganRandomParameter(organ_type):
        p.dxMin = 0.05
        p.dx = 0.1  # adjust resolution

# Simulate
plant.initialize()
for i in range(0, N):
    plant.simulate(dt)
    vp.write_plant("results/" + filename + "{:04d}".format(i), plant)
    print(i, '/', N)
