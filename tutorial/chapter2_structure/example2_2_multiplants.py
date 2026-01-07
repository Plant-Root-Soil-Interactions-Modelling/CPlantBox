"""multiple plants"""

import plantbox as pb
import plantbox.visualisation.vtk_plot as vp

path = path = "../../modelparameter/structural/plant/"
name = "fspm2023"

sim_time = 30  # days
n_colums = 3  # number of columns and rows #|\label{l2_2_2:Ncolsrows}|
dist = 40  # distance between the plants [cm] #|\label{l2_2_2:dist}|

# Initializes n_colums*n_colums plants
all_ = []
for i in range(0, n_colums):  # |\label{l2_2_2:iterationpositionbegin}|
    for j in range(0, n_colums):
        plant = pb.Plant()
        plant.readParameters(path + name + ".xml")
        seed = plant.getOrganRandomParameter(pb.seed)[0]
        seed.seedPos = pb.Vector3d(dist * i, dist * j, -3.0)  # cm
        plant.initialize(verbose = False)
        all_.append(plant)  # |\label{l2_2_2:iterationpositionend}|

# Simulate all plants
for plant in all_:  # |\label{l2_2_2:simbegin}|
    plant.simulate(sim_time, False)  # verbose = False #|\label{l2_2_2:simend}|

# Export results as single vtp files (as polylines)
ana = pb.SegmentAnalyser()  # |\label{l2_2_2:collectbegin}|
for i, plant in enumerate(all_):
    filename = "results/multplantsys_" + str(i)
    vp.write_plant(filename, plant)
    ana.addSegments(plant)  # collect all

# Write all into single file (as segments)
ana.write("results/multplantsys_all.vtp")  # |\label{l2_2_2:collectend}|
