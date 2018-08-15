import py_plantbox as pb
from rb_tools import *

# Simulate a root system
name = "sympodial_dichasium"
plant = pb.Plant()
#sca = pb.Organ.getScalar('organtype')
plant.openXML(name)






plant.initialize()
ana = pb.SegmentAnalyser(plant)

#for i in range (0,10):
#    plant.simulate(i)
plant.simulate(160) 


plant.write("../results/tree_sympodial_dichasium.vtp",15)
