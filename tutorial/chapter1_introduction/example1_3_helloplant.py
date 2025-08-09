import os

import plantbox as pb  # |\label{l13:cplantbox}|

# Optional plotting: works only if the repo's visualisation package and VTK are available
try:  # |\label{l13:vtk_plot}|
    import visualisation.vtk_plot as vp  # type: ignore

    _HAVE_VTK = True
except Exception:
    _HAVE_VTK = False

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

if _HAVE_VTK:
    # Interactive plot, using vtk, press x, y, z to change view, r to reset view, g to save png
    vp.plot_plant(plant, "age")  # e.g. organType, subType, age |\label{l13:plot_plant}|
