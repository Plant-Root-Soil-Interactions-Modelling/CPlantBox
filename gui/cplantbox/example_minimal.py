""" something basic"""
import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import visualisation.vtk_plot as vp
import numpy as np


def write_leafOnly():
    """ leafs on a  stem """
    plant = pb.Plant()

    seed = pb.SeedRandomParameter(plant)
    tap = pb.RootRandomParameter(plant)
    lateral = pb.RootRandomParameter(plant)  
    stem = pb.StemRandomParameter(plant)
    leaf = pb.LeafRandomParameter(plant)

    tap.name, lateral.name, stem.name, leaf.name = "Tap", "Lateral", "Stem", "Leaf"

    tap.subType = 1
    tap.lmax = 5
    tap.la, tap.lb, tap.ln = 1., 1., 1.
    tap.theta = 0.
    tap.r = 0.2 # 
    lateral.subType = 2
    lateral.lmax = 3
    lateral.la, lateral.lb, lateral.ln = 1, 1, 1
    lateral.theta = 90. / 180 * np.pi
    lateral.r = 0.1
    tap.successorOT = [[pb.root]]
    tap.successorST = [[2]]
    tap.successorP = [[1]]
    
    stem.subType = 1
    stem.lmax = 20
    stem.la, stem.lb, stem.ln = 3, 2, 3
    stem.theta = 0.
    stem.r = 1.
    stem.tropismT = 4
    stem.successorOT = [[pb.leaf]]
    stem.successorST = [[1]]
    stem.successorP = [[1]]

    leaf.parametrisationType = 0
    leaf.a = 0.02
    leaf.subType = 1
    leaf.lb = 1  # length of leaf stem
    leaf.la, leaf.lmax = 5, 11
    l_ = (leaf.lmax - leaf.lb) / 2  # radius of a round leaf
    leaf.areaMax = 0.5 * 3.145 * (l_ ** 2)
    leaf.leafGeometryPhi = np.array([-90, -45, 0., 45, 67.5, 70, 90]) / 180. * np.pi
    leaf.leafGeometryX = l_ * np.ones((7,))
    leaf.createGeometry()
    leaf.tropismS = 0.1

    plant.setOrganRandomParameter(seed)
    plant.setOrganRandomParameter(tap)
    plant.setOrganRandomParameter(lateral)
    plant.setOrganRandomParameter(stem)
    plant.setOrganRandomParameter(leaf)

    plant.initialize()
    plant.simulate(30)
    vp.plot_plant(plant, "organType")

    plant.writeParameters("params/leaf_only.xml")


def write_rootOnly():
    """ minimal root system """
    plant = pb.Plant()

    seed = pb.SeedRandomParameter(plant)
    tap = pb.RootRandomParameter(plant)
    lateral1 = pb.RootRandomParameter(plant)
    lateral2 = pb.RootRandomParameter(plant)
    basal = pb.RootRandomParameter(plant)
    stem = pb.StemRandomParameter(plant)

    tap.name, lateral1.name, lateral2.name, basal.name, stem.name = "Tap", "Lateral1", "Lateral2", "Basal", "Stem"

    seed.nC = 7

    tap.subType = 1
    tap.lmax = 60
    tap.la, tap.lb, tap.ln = 10, 2, 2
    tap.theta = 0.
    tap.tropismS = 0.2
    tap.tropismN = 2
    tap.r = 1.5
    tap.a = 0.15
    tap.successorOT = [[pb.root]]
    tap.successorST = [[2]]
    tap.successorP = [[1]]    
    lateral1.subType = 2
    lateral1.lmax = 25
    lateral1.la, lateral1.lb, lateral1.ln = 5, 1, 1
    lateral1.theta = 80. / 180 * np.pi
    lateral1.r = 1.
    lateral1.a = 0.1
    lateral1.successorOT = [[pb.root]]
    lateral1.successorST = [[3]]
    lateral1.successorP = [[1]]    
    lateral2.subType = 3
    lateral2.lmax = 10
    lateral2.la, lateral1.lb, lateral1.ln = 5, 1, 1
    lateral2.theta = 90. / 180 * np.pi
    lateral2.r = 0.2
    lateral2.a = 0.05
    
    basal.subType = 4
    basal.lmax = 20
    basal.theta = 70. / 180 * np.pi
    basal.la, basal.lb, basal.ln = 3, 2, 3



    stem.subType = 1
    stem.lmax = 3
    stem.la, stem.lb, stem.ln = 0., 1., 0.
    stem.theta = 0.
    stem.r = 1.
    stem.tropismT = 2

    plant.setOrganRandomParameter(seed)
    plant.setOrganRandomParameter(tap)
    plant.setOrganRandomParameter(lateral1)
    plant.setOrganRandomParameter(lateral2)
    plant.setOrganRandomParameter(basal)
    plant.setOrganRandomParameter(stem)

    plant.initialize()
    plant.simulate(50)
    vp.plot_plant(plant, "organType")

    plant.writeParameters("params/root_only.xml")


def write_stemOnly():
    """ minimal branched stem """
    plant = pb.Plant()

    seed = pb.SeedRandomParameter(plant)
    tap = pb.RootRandomParameter(plant)
    stem = pb.StemRandomParameter(plant)
    lateral = pb.RootRandomParameter(plant)
    lateralstem = pb.StemRandomParameter(plant)

    tap.name, lateral.name, stem.name, lateralstem.name = "Tap", "Lateral", "Main stem", "Lateral stem"

    tap.subType = 1
    tap.lmax = 5
    tap.la, tap.lb, tap.ln = 1., 1., 1.
    tap.theta = 0.
    tap.r = 0.2 
    lateral.subType = 2
    lateral.lmax = 3
    lateral.la, lateral.lb, lateral.ln = 1, 1, 1
    lateral.theta = 90. / 180 * np.pi
    lateral.r = 0.1
    tap.successorOT = [[pb.root]]
    tap.successorST = [[2]]
    tap.successorP = [[1]]

    stem.subType = 1
    stem.lmax = 20
    stem.la, stem.lb, stem.ln = 3, 2, 3
    stem.theta = 0.
    stem.r = 1.
    stem.tropismT = 2
    stem.betaDev = 10
    stem.a = 0.2
    lateralstem.a = 0.1
    lateralstem.subType = 2
    lateralstem.lmax = 10
    lateralstem.la, lateralstem.lb, lateralstem.ln = 1, 1, 1
    lateralstem.theta = 90. / 180 * np.pi
    lateralstem.tropismT = 4
    stem.successorOT = [[pb.stem]]
    stem.successorST = [[2]]
    stem.successorP = [[1]]

    plant.setOrganRandomParameter(seed)
    plant.setOrganRandomParameter(tap)
    plant.setOrganRandomParameter(lateral)
    plant.setOrganRandomParameter(stem)
    plant.setOrganRandomParameter(lateralstem)

    plant.initialize()
    plant.simulate(30)
    vp.plot_plant(plant, "organType")

    plant.writeParameters("params/stem_only.xml")


def update_maize(file_):
    plant = pb.Plant()
    plant.readParameters("params/" + file_)
    for p in plant.getOrganRandomParameter(pb.leaf):
        print(p.subType)
        p.theta = 10./ 180 * np.pi
        p.lb = 0  # length of leaf stem
        p.la, p.lmax = 49.12433414, 49.12433414
        p.areaMax = 71.95670914  # cm2, area reached when length = lmax
        NLeaf = 100
        phi = np.array([-90, -80, -45, 0., 45, 90]) / 180. * np.pi
        l = np.array([49.12433414, 1 , 1, 0.3, 1, 49.12433414])  # distance from leaf center
        p.tropismT = 1
        p.tropismN = 5
        p.tropismS = 0.05
        p.createLeafRadialGeometry(phi, l, NLeaf)
    
    plant.initialize()
    plant.simulate(30)
    vp.plot_plant(plant, "organType")
    
    # plant.writeParameters("params/P0.xml" + file_)

# write_leafOnly()
# write_rootOnly()
# write_stemOnly()

update_maize("P0.xml")
