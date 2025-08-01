from pathlib import Path
import sys

# project root = CPlantBox
root = Path(__file__).resolve().parents[2]

# add the built extension and Python helpers
sys.path.insert(0, str(root / "build" / "Release"))
sys.path.insert(0, str(root / "src"))

import plantbox as pb
print("Imported plantbox:", pb)



import plantbox as pb  # |\label{l13:cplantbox}|
import visualisation.vtk_plot as vp  # |\label{l13:vtk_plot}|

plant = pb.Plant()  # Create a new plant |\label{l13:plant}|

# Open plant and root parameter from a file
path = "../../modelparameter/structural/plant/"
name = "fspm2023"
plant.readParameters(path + name + ".xml")  # |\label{l13:readparameters}|

plant.initialize()  # Initialize |\label{l13:initialize}|

simtime = 40  # days
plant.simulate(simtime)  # Simulate|\label{l13:simulate}|

# Export final result (as vtp)
plant.write("results/example_plant.vtp")  # using polylines |\label{l13:write_poly}|

ana = pb.SegmentAnalyser(plant)
ana.write("results/example_plant_segs.vtp")  # using segments |\label{l13:write_segs}|

# Interactive plot, using vtk, press x, y, z to change view, r to reset view, g to save png
vp.plot_plant(plant, "age")  # e.g. organType, subType, age |\label{l13:plot_plant}|
