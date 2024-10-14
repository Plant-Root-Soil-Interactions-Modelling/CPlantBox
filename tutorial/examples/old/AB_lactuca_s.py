"""plant example

author : Arnaud Bouvry
date : 2022-03-24
adapted from "example_1e_plantShape"

"""
import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import visualisation.vtk_plot as vp

import numpy as np
import matplotlib.pyplot as plt

# #parameters for example:
adaptSeed = False

adaptStem = True

adaptLeaf = True
leafRadial = True  # radial or not
leafRound = True
leafLong = False
leafMaple = False
leafFlowershaped = False

anim = False
zoomLeafShape = False
export = False
getImage = False

# #create plant:
plant = pb.Plant()
# Open plant and root parameter from a file
path = "../../modelparameter/structural/plant/"
name = "Lactuca_s"
plant.readParameters(path + name + ".xml")

# #create rootsystem:
rs = pb.RootSystem()
# Open plant and root parameter from a file
path = "../../../modelparameter/plant/"
name = "Lactuca_s"
rs.readParameters(path + name + ".xml")

if adaptSeed:
    srp = pb.SeedRandomParameter(plant)  # with default values
    srp.firstTil = 0  # [day] first emergence of a tiler
    srp.delayTil = 0  # [day] delay between the emergence of tilers
    srp.maxTil = 0  # [-] number of tillers
    plant.setOrganRandomParameter(srp)

if adaptStem:
    for p in plant.getOrganRandomParameter(pb.stem):
        if (p.subType > 0):  # can be changed according to the suptypes of the plant
            p.nodalGrowth = 1  # < whether to implement the internodal growth
            p.delayLat = 1  # < delay between stem creation and start of nodal growth [day]
            # p.delayNG = 10  # < delay between lateral creation and growth [day]
            # p.tropismAge = 10 #< only used if tropsimT = 6
            plant.setOrganRandomParameter(p)

if adaptLeaf:
    for p in plant.getOrganRandomParameter(pb.leaf):

                # p.lmax - p.la - p.lb = leafMid = center of radial circle
        if (p.subType >= 2):  # leaf subtypes start at 2
            p.lb = 1  # length of leaf stem
            p.la, p.lmax = 3.5, 8.5
            p.areaMax = 10  # cm2, area reached when length = lmax
            N = 100  # N is rather high for testing

            if leafLong:
                # LongLeaf:
                p.lb = 1  # length of leaf stem
                p.la, p.lmax = 3.5, 8.5
                p.areaMax = 10  # cm2, area reached when length = lmax
                N = 100  # N is rather high for testing
                phi = np.array([-90, -45, 0., 45, 90]) / 180. * np.pi
                l = np.array([3, 2.2, 1.7, 2, 3.5])  # distance from leaf center

            if leafMaple:
                # Maple leaf:
                p.lb = 1  # length of leaf stem
                N = 100  # N is rather high for testing
                p.areaMax = 50
                p.la, p.lmax = 5, 11
                phi = np.array([-90, -45, 0., 45, 67.5, 70, 90]) / 180. * np.pi
                l = np.array([2, 2, 2, 4, 1, 1, 4])  # distance from leaf center

            if leafRound:
                # Round leaf:
                p.lb = 1  # length of leaf stem
                N = 100  # N is rather high for testing
                p.la, p.lmax = 5, 11
                # p.areaMax = 3.145*(((p.lmax-p.la - p.lb)/2)**2)
                # phi = np.array([-90, -45, 0., 45,67.5,70, 90]) / 180. * np.pi
                phi = np.arange(start = 0, stop = 2 * np.pi, step = np.pi / 4)
                l_ = (p.lmax - p.lb - p.la) / 2
                l = np.array([l_ for x_i in range(len(phi))])  # ([2, 2, 2, 4,1,1, 4]) #distance from leaf center

            if leafFlowershaped:
                # flower-shaped leaf:
                p.lb = 1  # length of leaf stem
                N = 100  # N is rather high for testing
                p.areaMax = 100
                p.la, p.lb, p.lmax = 5, 1, 11
                phi = np.array([-90., -67.5, -45, -22.5, 0, 22.5, 45, 67.5, 90]) / 180. * np.pi
                l = np.array([5., 1, 5, 1, 5, 1, 5, 1, 5])

            if leafRadial:
                p.createLeafRadialGeometry(phi, l, N)

            else:
                p.lb = 2  # length of leaf stem
                p.la, p.lmax = 3.5, 8.5
                p.areaMax = 10  # cm2, area reached when length = lmax
                N = 100  # N is rather high for testing
                y = np.array([-3, -3 * 0.7, 0., 3.5 * 0.7, 3.5])
                l = np.array([0., 2.2 * 0.7, 1.7, 1.8 * 0.7, 0.])
                p.createLeafGeometry(y, l, N)

            # p.tropismT = 4 # 6: Anti-gravitropism to gravitropism
            # p.tropismN = 1
            # p.tropismS = 0.1

            # p.tropismT = 6 # 6: Anti-gravitropism to gravitropism
            # p.tropismAge = 10 #< age at which tropism switch occures, only used if p.tropismT = 6
            # #p.tropismN = 5
            # #p.tropismS = 0.1

            plant.setOrganRandomParameter(p)

# snippet to display the parameters of various organ types
# for root subtypes
for p in plant.getOrganRandomParameter(pb.root):
    if (p.subType >= 0):
        print(p)

# for stem subtypes
for p in plant.getOrganRandomParameter(pb.stem):
    if (p.subType >= 0):
        print(p)

# for leaf
for p in plant.getOrganRandomParameter(pb.leaf):
    if (p.subType >= 2):
        print(p)

plant.initialize()
# rs.initialize()

plant.simulate(dt = 1, verbose = True)
plant.simulate(dt = 39, verbose = True)

# rs.simulate(dt=40, verbose=True)

vp.plot_plant(plant, "organType")
# vp.plot_plant(plant, "creationTime")

if getImage:
    # Simulate
    if not anim:
        plant.simulate(30, True)
    # Plot, using vtk
    vp.plot_plant(plant, "organType")
    # zoom on leaf--theory--2D
    print("2D leaf shape of a full grown leaf")
    lorg = plant.getOrgans(pb.leaf)[1]
    lrp = lorg.getLeafRandomParameter()
    leafRadial = (lrp.parametrisationType == 0)
    if leafRadial:
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

    # zoom on leaf--realized
    print("3D leaf shape of actual leaf")
    vp.plot_leaf(lorg)

