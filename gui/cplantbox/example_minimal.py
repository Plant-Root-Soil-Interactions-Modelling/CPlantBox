""" something basic"""
import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import visualisation.vtk_plot as vp
import numpy as np


def write_rootOnly():
    """ minimal root system """
    plant = pb.Plant()

    seed = pb.SeedRandomParameter(plant)
    tap = pb.RootRandomParameter(plant)
    lateral = pb.RootRandomParameter(plant)
    basal = pb.RootRandomParameter(plant)
    shootborne = pb.RootRandomParameter(plant)
    stem = pb.StemRandomParameter(plant)

    tap.name, lateral.name, basal.name, shootborne.name, stem.name = "Tap", "Lateral", "Basal", "Shootborne", "Stem"

    seed.nC = 7

    tap.subType = 1
    tap.lmax = 20
    tap.la, tap.lb, tap.ln = 3, 2, 3
    tap.theta = 0.
    tap.tropismS = 0.
    lateral.subType = 2
    lateral.lmax = 10
    lateral.la, lateral.lb, lateral.ln = 1, 1, 1
    lateral.theta = 90. / 180 * np.pi
    basal.subType = 4
    basal.lmax = 20
    basal.theta = 70. / 180 * np.pi
    basal.la, basal.lb, basal.ln = 3, 2, 3
    shootborne.subType = 5
    shootborne.lmax = 20
    shootborne.la, shootborne.lb, shootborne.ln = 3, 2, 3
    shootborne.theta = 80. / 180 * np.pi
    tap.successorOT = [[pb.root]]
    tap.successorST = [[2]]
    tap.successorP = [[1]]

    stem.subType = 1
    stem.lmax = 3
    stem.la, stem.lb, stem.ln = 0., 1., 0.
    stem.theta = 0.
    stem.r = 1.
    stem.tropismT = 2

    plant.setOrganRandomParameter(seed)
    plant.setOrganRandomParameter(tap)
    plant.setOrganRandomParameter(lateral)
    plant.setOrganRandomParameter(basal)
    plant.setOrganRandomParameter(shootborne)
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
    lateralstem = pb.StemRandomParameter(plant)

    tap.name, stem.name, lateralstem.name = "Tap", "Main stem", "Lateral stem"

    tap.subType = 1
    tap.lmax = 3
    tap.la, tap.lb, tap.ln = 0., 1., 0.
    tap.theta = 0.

    stem.subType = 1
    stem.lmax = 20
    stem.la, stem.lb, stem.ln = 3, 2, 3
    stem.theta = 0.
    stem.r = 1.
    stem.tropismT = 2
    stem.betaDev = 10
    stem.a = 1.
    lateralstem.a = 1
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
    plant.setOrganRandomParameter(stem)
    plant.setOrganRandomParameter(lateralstem)

    plant.initialize()
    plant.simulate(50)
    vp.plot_plant(plant, "organType")

    plant.writeParameters("params/stem_only.xml")


def write_leafOnly():
    """ leafs on a  stem """
    plant = pb.Plant()

    seed = pb.SeedRandomParameter(plant)
    tap = pb.RootRandomParameter(plant)
    stem = pb.StemRandomParameter(plant)
    leaf = pb.LeafRandomParameter(plant)

    tap.name, stem.name, leaf.name = "Tap", "Stem", "Leaf"

    tap.subType = 1
    tap.lmax = 3
    tap.la, tap.lb, tap.ln = 0., 1., 0.
    tap.theta = 0.

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
    leaf.areaMax = 3.145 * (l_ ** 2)
    leaf.leafGeometryPhi = np.array([-90, -45, 0., 45, 67.5, 70, 90]) / 180. * np.pi
    leaf.leafGeometryX = l_ * np.ones((7,))
    leaf.createGeometry()

    plant.setOrganRandomParameter(seed)
    plant.setOrganRandomParameter(tap)
    plant.setOrganRandomParameter(stem)
    plant.setOrganRandomParameter(leaf)

    plant.initialize()
    plant.simulate(50)
    vp.plot_plant(plant, "organType")

    plant.writeParameters("params/leaf_only.xml")


# write_rootOnly()
write_stemOnly()
# write_leafOnly()
