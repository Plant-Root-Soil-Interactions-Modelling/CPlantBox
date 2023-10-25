"""multiple root systems"""
import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import visualisation.vtk_plot as vp

path = path = "../../modelparameter/structural/plant/"
name = "fspm2023"

simtime = 90  # days
N = 3  # number of columns and rows
dist = 40  # distance between the root systems [cm]

# Initializes N*N root systems
all = []
for i in range(0, N):
    for j in range(0, N):
        plant = pb.Plant()
        plant.readParameters(path + name + ".xml")
        seed = plant.getOrganRandomParameter(pb.seed)[0]
        seed.seedPos = pb.Vector3d(dist * i, dist * j, -3.)  # cm
        plant.initialize(verbose = False)
        all.append(plant)

# Simulate all plants
for plant in all:
    plant.simulate(simtime, False)  # verbose = False

# Export results as single vtp files (as polylines)
ana = pb.SegmentAnalyser()
for i, plant in enumerate(all):
      vtpname = "results/topics_virtual4_" + str(i) + ".vtp"
      plant.write(vtpname)
      ana.addSegments(plant)  # collect all

# Write all into single file (as segments)
ana.write("results/topics_virtual4_all.vtp")

# Plot, using vtk
vp.plot_roots(ana, "radius")  # TODO plot_plants is not working for the SegmentAnalyser
