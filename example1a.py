import py_plantbox as pb
from rb_tools import *



plant = pb.Plant()

# Open  and root parameter from a file
name = "CPlantBox_test_leaf_tree" # "Anagallis_femina_Leitner_2010"
plant.openXML(name)

# Initialize
plant.initialize()

# Simulate
plant.simulate(60, True)

# Export final result (as vtp)
#2 = root
#4 = stem
#8 = leaf
#15 = all

plant.write("results/example_1a.vtp", 15)
