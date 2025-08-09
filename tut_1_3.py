import os

import plantbox as pb  # |\label{l13:cplantbox}|
import vtk_plot as vp  # type: ignore

plant = pb.Plant()  # Create a new plant |\label{l13:plant}|

# Open plant and root parameter from a file using packaged data paths
path = os.path.join(pb.data_path(), "structural", "plant") + "/"
name = "fspm2023"
plant.readParameters(path + name + ".xml")  # |\label{l13:readparameters}|

plant.initialize()  # Initialize |\label{l13:initialize}|

simtime = 40  # days
plant.simulate(simtime)  # Simulate|\label{l13:simulate}|

os.makedirs("results", exist_ok=True)
# Export final result (as vtp)
plant.write("results/example_plant.vtp")  # using polylines |\label{l13:write_poly}|

ana = pb.SegmentAnalyser(plant)
ana.write("results/example_plant_segs.vtp")  # using segments |\label{l13:write_segs}|

vp.plot_plant(plant, "age")  # e.g. organType, subType, age |\label{l13:plot_plant}|
