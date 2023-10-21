
import sys; sys.path.append("../.."); sys.path.append("../../src/")
path =  "../../modelparameter/structural/plant/"
sys.path.append( path)
import plantbox as pb
import visualisation.vtk_plot as vp
from example1f import template_text

p = pb.MappedPlant(2)
import os
print(os.path.isfile("P0_plant_new.xml"))
p.readParameters("P0_plant_new.xml") #the parameters are given directly as string and not via a text file        
p.initialize(False)
time = 100
p.simulate(time, False)
ana = pb.SegmentAnalyser()  # see example 3b

#vp.plot_plant(p, "organType")
p.write("testBranching.vtp")

