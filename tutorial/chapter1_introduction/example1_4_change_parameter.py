import sys;sys.path.append("../..");sys.path.append("../../src/")  # |\label{l14:path}|

import plantbox as pb  # |\label{l14:cplantbox}|
import visualisation.vtk_plot as vp  # |\label{l14:vtk_plot}|

plant = pb.Plant()  # Create a new plant |\label{l14:plant}|

# Open plant and root parameter from a file
path = "../../modelparameter/structural/plant/"
name = "fspm2023"
plant.readParameters(path + name + ".xml")  # |\label{l14:readparameters}|

# Get input parameter by organ type
root = plant.getOrganRandomParameter(pb.root)  # |\label{l14:root}|
stem = plant.getOrganRandomParameter(pb.stem)  # |\label{l14:stem}|
leaf = plant.getOrganRandomParameter(pb.leaf)  # |\label{l14:leaf}|
seed = plant.getOrganRandomParameter(pb.seed)  # |\label{l14:seed}|

print(root, stem, leaf, seed)  # Print seed, root, stem, leaf, and seed parameters |\label{l14:print}|

# Change a parameter
root[1].r = 5  # Change elongation rate (r [cm/day]) |\label{l14:change_params_r}|
root[1].ln = 0.25  # Change inter-lateral distance (ln [cm]) |\label{l14:change_params_ln}|
root[2].dx = 0.5  # Change axial resolution (dx [cm]) |\label{l14:change_params_dx}|
root[2].dxMin = 0.1  # Change minimal axial resolution (dxMin [cm])

print(root[1])  # Print new root parameters |\label{l14:print_new}|

plant.initialize()  # Initialize |\label{l14:initialize}|

simtime = 40  # days
plant.simulate(simtime)  # Simulate|\label{l14:simulate}|

# Plot
vp.plot_plant(plant, "organType")
