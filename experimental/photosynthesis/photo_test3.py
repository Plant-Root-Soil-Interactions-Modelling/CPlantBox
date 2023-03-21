import sys; sys.path.append("../.."); sys.path.append("../../src")
import plantbox as pb
from functional.xylem_flux import XylemFluxPython  # Python hybrid solver
import numpy as np
import matplotlib.pyplot as plt
from functional.photosynthesis_cpp import PhotosynthesisPython
import visualisation.vtk_plot as vp # for quick 3d vizualisations

simtime = 14  # [day] 
pl = pb.MappedPlant() #for plant objects
path = "../../modelparameter/structural/plant/" 
name = "fspm2023" 
pl.readParameters(path + name + ".xml")

soilSpace = pb.SDF_PlantContainer(500,500, 500, True) #to avoid root growing aboveground
pl.setGeometry(soilSpace) # creates soil space to stop roots from growing out of the soil

pl.initialize()
pl.simulate(simtime)
vp.plot_plant(pl, "organType")