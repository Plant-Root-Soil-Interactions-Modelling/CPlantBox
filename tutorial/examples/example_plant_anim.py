"""plant example"""
import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import visualisation.vtk_plot as vp

import numpy as np

plant = pb.Plant()

# Open plant and root parameter from a file
path = "../../modelparameter/structural/plant/"
name = "0"  # CPlantBox_test_leaf_tree00

# LEAFS smallPlant_mgiraud "manyleaves"
# NO LEAFS "CPlantBox_test_leaf_tree22"  # "chicon_entire"  # "Anagallis_femina_leaf_shape"  # "Anagallis_femina_Leitner_2010"

# BREAKS MY COMPUTER Swiss_chard

plant.readParameters(path + name + ".xml")


# soil = pb.SDF_PlantContainer(1.e6, 1.e6, 1.e6, False)
# plant.setGeometry(soil)

# increase resolution
for p in plant.getOrganRandomParameter(pb.root):
    p.dx = 0.2

# Initialize
plant.initialize()

dt = 1
N = 400
min_ = np.array([-20, -20, -50])
max_ = np.array([20, 20, 30.])

# test = plant.getOrgans(pb.leaf)
# print("test")

anim = vp.AnimateRoots(plant)
anim.min = min_
anim.max = max_
anim.res = [1, 1, 1]
anim.file = "results/example_plant"
anim.avi_name = "results/example_"
anim.plant = True
anim.start()

for i in range(0, N):

    plant.simulate(dt, False)
    anim.root_name = "organType"
    anim.update()

