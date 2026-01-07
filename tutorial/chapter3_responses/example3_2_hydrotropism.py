"""hydrotropism in a thin layer"""

import plantbox as pb
import plantbox.visualisation.vtk_plot as vp

plant = pb.Plant()
path = "../../modelparameter/structural/rootsystem/"  # |\label{l3_2_hydrotropism:libsstart}|
filename = "Anagallis_femina_Leitner_2010"
plant.readParameters(path + filename + ".xml")  # |\label{l3_2_hydrotropism:libsend}|

# Manually set tropism to hydrotropism for the first ten root types
sigma = [0.4, 1.0, 1.0, 1.0, 1.0] * 2  # |\label{l3_2_hydrotropism:tsetstart}|
for p in plant.getOrganRandomParameter(pb.root):
    p.dx = 0.25  # adjust resolution
    p.tropismT = pb.TropismType.hydro
    p.tropismN = 1.8  # strength of tropism
    p.tropismS = sigma[p.subType - 1]  # |\label{l3_2_hydrotropism:tsetend}|

# Static soil property in a thin layer
maxS = 0.7  # maximal # |\label{l3_2_hydrotropism:soilpropstart}|
minS = 0.1  # minimal
slope = 5  # linear gradient between min and max (cm)
box = pb.SDF_PlantBox(30, 30, 2)  # cm
layer = pb.SDF_RotateTranslate(box, pb.Vector3d(0, 0, -16))
soil_prop = pb.SoilLookUpSDF(layer, maxS, minS, slope)  # |\label{l3_2_hydrotropism:soilpropend}|

# Set the soil properties before calling initialize
plant.setSoil(soil_prop)  # |\label{l3_2_hydrotropism:set&simstart}|

# Initialize
plant.initialize()

# Simulate
sim_time = 100  # e.g. 30 or 60 days
dt = 1
n_steps = round(sim_time / dt)
for _ in range(0, n_steps):
    # in a dynamic soil setting you would need to update the soil properties (soil_prop)
    plant.simulate(dt)  # |\label{l3_2_hydrotropism:set&simend}|

# Export results (as vtp) # |\label{l3_2_hydrotropism:resultsstart}|
plant.write("results/example3_2_hydrotropism.vtp")

# Export geometry of static soil
plant.setGeometry(layer)  # just for vizualisation
plant.write("results/example3_2_hydrotropism.py")

# Plot, using vtk
vp.plot_roots_and_container(plant, layer, "type")  # |\label{l3_2_hydrotropism:resultsend}|
