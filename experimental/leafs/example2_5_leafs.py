import numpy as np

import plantbox as pb
import plantbox.visualisation.vtk_plot as vp


def mint(leaf_rp):
    """Create a plant with a mint leaf shape"""
    leaf_rp.parametrisationType = 0  # 2D shape type : 0 .. radial, 1..along main axis
    leaf_rp.lb = 1.5  # petiole length
    leaf_rp.la = 3.5  # mid point for radial disretization is at 3.5 cm from tip
    leaf_rp.lmax = 7.5  #  maximal length including petiole
    leaf_rp.areaMax = 20  # [cm2] maximal leaf area, blade width results from area and length
    phi = np.array([-90, -45, 0.0, 45, 90]) / 180.0 * np.pi
    l = np.array([3, 2.2, 1.7, 2, 3.5])
    leaf_rp.createLeafRadialGeometry(phi, l, 100)


def parsley(leaf_rp):
    """Create a plant with a pasley leaf shape"""
    leaf_rp.parametrisationType = 0  # 2D shape type : 0 .. radial, 1..along main axis
    leaf_rp.lb = 1  # petiole length
    leaf_rp.la = 5.4  # mid point for radial disretization is at 5.4 cm from tip
    leaf_rp.lmax = 5.4 + leaf_rp.lb + 0.8  # maximal length including petiole
    leaf_rp.areaMax = 40  # [cm2] maximal leaf area, blade width results from area and length
    phi = np.array([-90, -45, -25, -15, 0, 10, 20, 30, 40, 45, 55, 60, 65, 90]) / 180.0 * np.pi
    l = np.array([0.8, 0.5, 1.5, 2, 2.6, 2, 2.6, 2.7, 3.1, 3.3, 2.6, 3.3, 3.8, 5.4])
    leaf_rp.createLeafRadialGeometry(phi, l, 100)


def example(leaf_rp):
    """Create a plant with a star shaped leaf shape"""
    leaf_rp.parametrisationType = 0  # 2D shape type : 0 .. radial, 1..along main axis
    leaf_rp.lb = 1.0  # petiole length
    leaf_rp.la = 5.0  # mid point for radial disretization is at 5.0 cm from tip
    leaf_rp.lmax = 11.0  # maximal length including petiole
    leaf_rp.areaMax = 50  # [cm2] maximal leaf area, blade width results from area and length
    phi = np.array([-90.0, -67.5, -45, -22.5, 0, 22.5, 45, 67.5, 90]) / 180.0 * np.pi
    l = np.array([5.0, 2, 5, 2, 5, 2, 5, 2, 5])
    leaf_rp.createLeafRadialGeometry(phi, l, 100)


def long_leaf(leaf_rp):
    """Create a plant with a long leaf shape"""
    leaf_rp.parametrisationType = 1  # 2D shape type : 0 .. radial, 1..along main axis
    leaf_rp.areaMax = 10  # cm2, area reached when length = lmax, blade width results from area and length
    leaf_rp.lb = 1  # cm
    y = np.array([-3, -3 * 0.7, 0.0, 3.5 * 0.7, 3.5])
    leaf_rp.la = y[-1]  # leaf blade mid is located lmax-la
    leaf_rp.lmax = leaf_rp.lb + (leaf_rp.la - y[0])  # petiole length + blade lenght
    l = np.array([0.0, 2.2 * 0.7, 1.7, 1.8 * 0.7, 0.0])
    leaf_rp.createLeafGeometry(y, l, 100)


if __name__ == "__main__":

    plant = pb.Plant()
    seed_rp = pb.SeedRandomParameter(plant)
    seed_rp.delayDefinitionShoot = pb.dd_distance
    plant.setOrganRandomParameter(seed_rp)

    root_rp = pb.RootRandomParameter(plant)
    root_rp.subType = 1
    # root_rp.lmax = 1
    # root_rp.r = 0.1
    # root_rp.theta = 0
    plant.setOrganRandomParameter(root_rp)

    stem_rp = pb.StemRandomParameter(plant)
    stem_rp.subType = 1
    stem_rp.la = 2.0
    stem_rp.ln = 2.0
    stem_rp.lb = 2
    stem_rp.lmax = 6  # 7.5  # 8.0
    stem_rp.r = 0.2
    stem_rp.successor = [[2]]
    stem_rp.successorP = [[1]]
    stem_rp.successorOT = [[pb.leaf]]
    stem_rp.successorNo = [1]
    stem_rp.ldelay = 4

    print(stem_rp.dx, stem_rp.dxMin)
    plant.setOrganRandomParameter(stem_rp)

    leaf_rp = pb.LeafRandomParameter(plant)
    leaf_rp.shapeType = 2  # Shape of the leaf: 0: cylinder (a = radius), 1: cuboid (a = thickness, Width_blade, Width_petiole), 2: 2D (leafGeometryPhi, leafGeometryX, areaMax)
    leaf_rp.subType = 2
    leaf_rp.tropismS = 0.0
    leaf_rp.theta = 70.0 / (2 * np.pi)
    leaf_rp.dx = 0.1
    leaf_rp.a = 0.05
    leaf_rp.rotBeta = np.pi  #####################
    leaf_rp.betaDev = 0.0
    leaf_rp.initBeta = 0.0
    leaf_rp.ldelay = 4  # parent decides delay for leaf laterals
    leaf_rp.r = 1  # growth rate of centerline [cm/day]

    # parsley(leaf_rp)
    long_leaf(leaf_rp)

    plant.setOrganRandomParameter(leaf_rp)

    plant.initialize()
    plant.simulate(22)
    vp.plot_plant(plant, "subType")  # creationTime


# problems:
# Two leafes, e.g.  when lmax is 7.5
# rotBeta not working (?)
