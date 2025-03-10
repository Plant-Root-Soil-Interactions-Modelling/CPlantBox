""" something on PlantVisualiser, cpbvis """
import sys; sys.path.append("../.."); sys.path.append("../../src/"); sys.path.append("./"); sys.path.append("./src/")

import plantbox as pb
import visualisation.vtk_plot as vp
import visualisation.vis_tools as cpbvis

import numpy as np

filename = "../../modelparameter/structural/plant/shootOnly_fspm2023.xml"
 

time_init = 10
time_max = 20
dt = 1
N = int((time_max -time_init  )/dt)

leaf_res = 30
# create a plant
plant = pb.MappedPlant()
plant.readParameters(filename)
#vis = pb.PlantVisualiser(plant)
plant.initialize()
#vis.SetGeometryResolution(8)
#vis.SetLeafResolution(leaf_res)
plant.simulate(time_init, False)

for i in range(N):
    print('loop',i+1,'/',N)
    plant.simulate(dt, False)
    

    #vis.ResetGeometry()
    #vis.ComputeGeometryForOrganType(pb.stem, False)
    #vis.ComputeGeometryForOrganType(pb.leaf, False)
    
    # Write the geometry to file#
    #data = cpbvis.PolydataFromPlantGeometry(vis)
    #cpbvis.WritePolydataToFile(data, "./results/vis_plant_option1_at"+str(i) ) # from vis
    vp.write_plant("./results/vis_plant_"+str(i), plant)  # will write XXX.vtp (with centerlines), and XXX_leafs.vtp (as polygones)

vp.plot_plant(plant, "subType")

print("fin")


'''

'''