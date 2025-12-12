"""multiple root systems"""

import plantbox as pb
import plantbox.visualisation.vtk_plot as vp

path = path = "../../modelparameter/structural/plant/"
name = "fspm2023"

simtime = 30  # days
N = 3  # number of columns and rows #|\label{l2_2_2:Ncolsrows}|
dist = 40  # distance between the root systems [cm] #|\label{l2_2_2:dist}|

# Initializes N*N root systems
all_ = []
for i in range(0, N):  # |\label{l2_2_2:iterationpositionbegin}|
    for j in range(0, N):
        plant = pb.Plant()
        plant.readParameters(path + name + ".xml")
        seed = plant.getOrganRandomParameter(pb.seed)[0]
        seed.seedPos = pb.Vector3d(dist * i, dist * j, -3.0)  # cm
        plant.initialize(verbose=False)
        all_.append(plant)  # |\label{l2_2_2:iterationpositionend}|

# Simulate all plants
for plant in all_:  # |\label{l2_2_2:simbegin}|
    plant.simulate(simtime, False)  # verbose = False #|\label{l2_2_2:simend}|

# Export results as single vtp files (as polylines)
ana = pb.SegmentAnalyser()  # |\label{l2_2_2:collectbegin}|
for i, plant in enumerate(all_):
    filename = "results/multrootsys_" + str(i)
    vp.write_plant(filename, plant)
    ana.addSegments(plant)  # collect all

# Write all into single file (as segments)
ana.write("results/multrootsys_all.vtp")  # |\label{l2_2_2:collectend}|
