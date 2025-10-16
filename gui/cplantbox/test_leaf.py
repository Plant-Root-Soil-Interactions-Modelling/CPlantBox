"""plant example

adapted from "example_plant_anim", "example10_nodalGrowth"
"CPlantBox_tutorial_FSPM2020.ipynb", "test_leafparameter"
We present here parameters specific for leaves and stems

"""
import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import visualisation.vtk_plot as vp

import numpy as np
import matplotlib.pyplot as plt

adaptLeaf = False
getImage = True

plant = pb.Plant()

# # Open plant and root parameter from a file
# path = "../../modelparameter/structural/plant/"
# name = "example1e"
# plant.readParameters(path + name + ".xml")

path = "params/"
# plant.readParameters(path +  "leaf_only.xml")
# plant.readParameters(path +  "root_only.xml")
# plant.readParameters(path +  "stem_only.xml")
# plant.readParameters(path +  "P0.xml")
# plant.readParameters(path +  "Triticum_aestivum_test_2021.xml")
plant.readParameters(path +  "fspm2023.xml")


# if adaptLeaf:
#     for p in plant.getOrganRandomParameter(pb.leaf):
#
#                 # p.lmax - p.la - p.lb = leafMid = center of radial circle
#         if (p.subType >= 2):  # leaf subtypes start at 2
#             p.lb = 1  # length of leaf stem
#             p.la, p.lmax = 3.5, 8.5
#             p.areaMax = 10  # cm2, area reached when length = lmax
#             N = 100  # N is rather high for testing
#             if leafRadial:
#
#                 # LongLeaf:
#                 p.lb = 1  # length of leaf stem
#                 p.la, p.lmax = 3.5, 8.5
#                 p.areaMax = 10  # cm2, area reached when length = lmax
#                 N = 100  # N is rather high for testing
#                 phi = np.array([-90, -45, 0., 45, 90]) / 180. * np.pi
#                 l = np.array([3, 2.2, 1.7, 2, 3.5])  # distance from leaf center
#
#                 # Maple leaf:
#                 p.lb = 1  # length of leaf stem
#                 N = 100  # N is rather high for testing
#                 p.areaMax = 50
#                 p.la, p.lmax = 5, 11
#                 phi = np.array([-90, -45, 0., 45, 67.5, 70, 90]) / 180. * np.pi
#                 l = np.array([2, 2, 2, 4, 1, 1, 4])  # distance from leaf center
#
#                 # Round leaf:
#                 p.lb = 1  # length of leaf stem
#                 N = 100  # N is rather high for testing
#                 p.la, p.lmax = 5, 11
#                 p.areaMax = 3.145 * (((p.lmax - p.la - p.lb) / 2) ** 2)
#                 phi = np.array([-90, -45, 0., 45, 67.5, 70, 90]) / 180. * np.pi
#                 l_ = (p.lmax - p.lb - p.la) / 2
#                 l = np.array([l_ for x_i in range(len(phi))])  # ([2, 2, 2, 4,1,1, 4]) #distance from leaf center
#
#                 # flower-shaped leaf:
#                 # p.lb =  1 # length of leaf stem
#                 # N = 100  # N is rather high for testing
#                 # p.areaMax = 100
#                 # p.la, p.lb, p.lmax = 5, 1, 11
#                 # phi = np.array([-90., -67.5, -45, -22.5, 0, 22.5, 45, 67.5, 90]) / 180. * np.pi
#                 # l = np.array([5., 1, 5, 1, 5, 1, 5, 1, 5])
#
#                 p.createLeafRadialGeometry(phi, l, N)
#
#             else:
#                 p.lb = 2  # length of leaf stem
#                 p.la, p.lmax = 3.5, 8.5
#                 p.areaMax = 10  # cm2, area reached when length = lmax
#                 N = 100  # N is rather high for testing
#                 y = np.array([-3, -3 * 0.7, 0., 3.5 * 0.7, 3.5])
#                 l = np.array([0., 2.2 * 0.7, 1.7, 1.8 * 0.7, 0.])
#                 p.createLeafGeometry(y, l, N)
#
#             p.tropismT = 6  # 6: Anti-gravitropism to gravitropism
#             # p.tropismN = 5
#             # p.tropismS = 0.1
#             p.tropismAge = 10  # < age at which tropism switch occures, only used if p.tropismT = 6
#             plant.setOrganRandomParameter(p)

for p in plant.getOrganRandomParameter(pb.leaf):
    if (p.subType >= 1):
        print("subType", p.subType)
        print("shapeType", p.shapeType)
        print("parametrisationType", p.parametrisationType)        
        print("leafGeometryPhi", p.leafGeometryPhi)
        print("leafGeometryX", p.leafGeometryX)
        print("areaMax", p.areaMax)
        print("la", p.la)
        print("lb", p.lb)


plant.initialize()

if getImage:
    # Simulate
    plant.simulate(20, True)
    # Plot, using vtk
    vp.plot_plant(plant, "organType")
    # zoom on leaf--theory--2D
    print("2D leaf shape of a full grown leaf")
    lorg = plant.getOrgans(pb.leaf)[1]
    lrp = lorg.getLeafRandomParameter()    
    leafRadial = (lrp.parametrisationType == 0)
    
    if leafRadial:
        print("radial")
        print(lrp.leafGeometry)
        N = len(lrp.leafGeometry)        
        yy = np.linspace(0, lorg.leafLength(), N)        
        geom_x, geom_y = [], []
        for i, x in enumerate(lrp.leafGeometry):
            geom_x.extend(x)
            geom_y.extend([yy[i]] * len(x))
        geom_x = np.array(geom_x)
        geom_y = np.array(geom_y)
        a = lorg.leafArea() / lorg.leafLength()  # scale radius
        plt.plot(geom_x * a, geom_y, "g*")
        plt.plot(-geom_x * a, geom_y, "g*")

    else:
        print("axial")
        print(lrp.leafGeometry)
        geom_x_a = np.array([0])
        geom_x_b = np.array([ x[-1] for x in lrp.leafGeometry])  # normalized x value along length
        geom_x = np.concatenate((geom_x_a, geom_x_b))
        geom_y_a = np.array([0])
        geom_y_b = np.linspace(lrp.lb, lorg.leafLength() + lrp.lb, len(geom_x_b))
        geom_y = np.concatenate((geom_y_a, geom_y_b))
        a = lorg.leafArea() / lorg.leafLength()  # scale radius
        plt.plot(geom_x * a, geom_y, "g-*")
        plt.plot(-geom_x * a, geom_y, "g-*")
        
    plt.ylim([0, lrp.lmax + 1])
    plt.xlim([-a - 1, a + 1])
    plt.axis('scaled')
    plt.show()

