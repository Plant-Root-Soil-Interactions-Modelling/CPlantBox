""" introctionary example """
import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import visualisation.vtk_plot as vp

import numpy as np

# Create a new plant
plant = pb.Plant()

# Open plant and root parameter from a file
path = "../../modelparameter/structural/plant/"
name = "hello_world"
plant.readParameters(path + name + ".xml")

# Initialize
plant.initialize()

# Simulate
plant.simulate(40, verbose = True)

# Export final result (as vtp)
plant.write("results/example_plant.vtp")  # using polylines

ana = pb.SegmentAnalyser(plant)
ana.write("results/example_plant_segs.vtp")  # using segments

# Interactive plot, using vtk, press x, y, z to change view, r to reset view, g to save png
vp.plot_plant(plant, "subType")
