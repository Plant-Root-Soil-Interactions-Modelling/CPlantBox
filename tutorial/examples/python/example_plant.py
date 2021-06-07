"""plant example"""
import sys
sys.path.append("../../.."); sys.path.append("../../../src/python_modules")
import plantbox as pb
import vtk_plot as vp

plant = pb.Plant()

# Open plant and root parameter from a file
path = "../../../modelparameter/plant/"
name = "Swiss_chard"   

# LEAFS smallPlant_mgiraud "manyleaves"
# NO LEAFS "CPlantBox_test_leaf_tree22"  # "chicon_entire"  # "Anagallis_femina_leaf_shape"  # "Anagallis_femina_Leitner_2010"

# BREAKS MY COMPUTER Swiss_chard

plant.readParameters(path + name + ".xml")

# print radii
print("leafs")
for p in plant.getOrganRandomParameter(pb.leaf):
    print(p.subType, p.a)
    p.a = 0.05
    p.a_s = 0
print("stem")
for p in plant.getOrganRandomParameter(pb.stem):
    print(p.subType, p.a)
    p.a = 0.2
    p.a_s = 0
print("roots")
for p in plant.getOrganRandomParameter(pb.root):
    print(p.subType, p.a)
    p.a = 0.05
    p.a_s = 0

soil = pb.SDF_PlantContainer(1.e6, 1.e6, 1.e6, False)
# plant.setGeometry(soil)

# increase resolution
for p in plant.getOrganRandomParameter(pb.root):
    p.dx = 0.2

# Initialize
plant.initialize()

# Simulate
plant.simulate(90, True)

# Export final result (as vtp)
plant.write("results/example_plant.vtp")

ana = pb.SegmentAnalyser(plant)
ana.write("results/example_plant_segs.vtp")

# Plot, using vtk
vp.plot_roots(plant, "organType")  # "creationTime"
