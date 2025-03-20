"""
growing plant example, use of AnimateRoots 

TODO we could add how to create avi form png (e.g. on linux), and remove files again
"""
import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import visualisation.vtk_plot as vp

import numpy as np

plant = pb.MappedPlant(0)

# Open plant and root parameter from a file
path = "../../modelparameter/structural/plant/"
name = "Triticum_aestivum_test_2021_4POF"  # CPlantBox_test_leaf_tree00

# LEAFS smallPlant_mgiraud "manyleaves"
# NO LEAFS "CPlantBox_test_leaf_tree22"  # "chicon_entire"  # "Anagallis_femina_leaf_shape"  # "Anagallis_femina_Leitner_2010"

# BREAKS MY COMPUTER Swiss_chard

plant.readParameters(path + name + ".xml")

for p in plant.getOrganRandomParameter(pb.leaf):
    p.lb = 2 # length of leaf stem
    p.la,  p.lmax = 49.12433414, 49.12433414
    p.areaMax = 71.95670914  # cm2, area reached when length = lmax
    NLeaf = 100  
    phi = np.array([-90,-80, -45, 0., 45, 90]) / 180. * np.pi
    l = np.array([49.12433414,1 ,1, 0.3, 1, 49.12433414]) #distance from leaf center
    #p.tropismT = 1 # 6: Anti-gravitropism to gravitropism
    #p.tropismN = 5
    #p.tropismS = 0.05
    
    p.createLeafRadialGeometry(phi, l, NLeaf)
# print ra
plant.initialize()

dt = 1
N = int(np.ceil(25/dt))
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
nodesOLd = np.array([])
nodesOLd2 = np.array([])
for i in range(0, N):

    plant.simulate(dt, False)
    if len(plant.getOrgans(4))>0:
        print('leaf0 id',plant.getOrgans(4)[0].getId())
    if False:
        print(plant.getSimTime(),'len(plant.getOrgans(4))',len(plant.getOrgans(4)))
        if len(plant.getOrgans(4))>0:
            leafs = plant.getOrgans(4)[0]
            nodesNew= np.array([[leafs.getNode(ii).x,leafs.getNode(ii).y,leafs.getNode(ii).z] for ii in range(leafs.getNumberOfNodes())])
            segLength =np.diff( np.array([leafs.getLength(ii) for ii in range(leafs.getNumberOfNodes())]))
            if(len(nodesOLd)>0):
                print('leaf1', nodesNew,nodesOLd - nodesNew[:len(nodesOLd)])
                print('segLength',segLength)
            nodesOLd = nodesNew
        for ii, leaf in enumerate(plant.getOrgans(4)):
            print(ii, leaf.getId(),np.array([[leaf.getNode(ii).x,leaf.getNode(ii).y,leaf.getNode(ii).z] for ii in range(leaf.getNumberOfNodes())]))
        if False:#len(plant.getOrgans(4))>1:
            leafs = plant.getOrgans(4)[1]
            nodesNew2= np.array([[leafs.getNode(ii).x,leafs.getNode(ii).y,leafs.getNode(ii).z] for ii in range(leafs.getNumberOfNodes())])
            if(len(nodesOLd2)>0):
                print('leaf2', nodesOLd2 - nodesNew2[:len(nodesOLd2)])
            nodesOLd2 = nodesNew2
    anim.root_name = "organType"
    anim.update()

