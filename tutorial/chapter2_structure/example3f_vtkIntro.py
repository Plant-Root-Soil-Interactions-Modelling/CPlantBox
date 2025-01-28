"""increase axial resolution (e.g. for animation)"""
import sys; sys.path.append("../.."); sys.path.append("../../src/")
import plantbox as pb
import visualisation.vtk_plot as vp
import numpy as np


""" plant """  # |\label{3f:plantStart}|
plant = pb.MappedPlant(0)
path = "../../modelparameter/structural/plant/"
name = "fspm2023" 
plant.readParameters(path + name + ".xml")
simtime = 60. # days

# Simulate
plant.initialize()
plant.simulate(simtime)  # |\label{3f:plantEnd}|

# plot results
vp.plot_plant(plant, "subType") # |\label{3f:option1}|
ana = pb.SegmentAnalyser(plant.mappedSegments())
vp.plot_plant(ana,"creationTime") # |\label{3f:option2}|

""" add parameters """ 
random_array = np.random.rand(len(plant.segments)) 
ana.addData("random_array",random_array) # |\label{3f:addData}|
vp.plot_plant(ana,"random_array")


""" root system """  # |\label{3f:rootsystem}|
plant = pb.MappedPlant()  
path = "../../modelparameter/structural/rootsystem/"
name = "Anagallis_femina_Leitner_2010"
plant.readParameters(path + name + ".xml")
plant.initialize()
plant.simulate(simtime)  # |\label{l41:rootsystem_end}|

""" periodic representation """ 
ana = pb.SegmentAnalyser(plant.mappedSegments())
ana.mapPeriodic(5, 5)           # |\label{3f:mapPeriodic}|
vp.plot_plant(ana, "creationTime")

