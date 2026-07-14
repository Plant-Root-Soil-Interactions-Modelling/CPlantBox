import numpy as np

import plantbox as pb
import plantbox.visualisation.vtk_plot as vp


class EmptyPlant(pb.Plant):
    """adds functionality to initialize with an empty plant, where organs are added later"""

    def __init__(self, seed_num=0):
        super().__init__(seed_num)

    def initialize_empty(self):
        """initializes an empt plant, use addOrgan to later add organs, e.g. leaves, roots, stems, etc.
        called instead of initialize(), initializeLB(), or initializeDB()"""
        self.setOrganRandomParameter(pb.SeedRandomParameter(self))
        self.reset()  # resets plant
        seed = pb.Seed(self)  # create a seed organ with the given random parameter
        self.addOrgan(seed)  # add the seed organ to the plant
        self.oldNumberOfNodes = self.getNumberOfNodes()
        self.initCallbacks()


plant = EmptyPlant()


leaf_rp = pb.LeafRandomParameter(plant)
leaf_rp.parametrisationType = 0  # 2D shape type : 0 .. radial, 1..along main axis
leaf_rp.shapeType = 2  # Shape of the leaf: 0: cylinder (a = radius), 1: cuboid (a = thickness, Width_blade, Width_petiole), 2: 2D (leafGeometryPhi, leafGeometryX, areaMax)
leaf_rp.subType = 1
leaf_rp.lmax = 10
leaf_rp.tropismS = 0.0
leaf_rp.theta = 90.0 / (2 * np.pi)

# Mint
leaf_rp.la, leaf_rp.lb, leaf_rp.lmax, leaf_rp.ln, leaf_rp.r, leaf_rp.dx = 3.5, 1.0, 7.5, 3, 1, 0.5
phi = np.array([-90, -45, 0.0, 45, 90]) / 180.0 * np.pi
l = np.array([3, 2.2, 1.7, 2, 3.5])
N = 30  # N is rather high for testing
leaf_rp.createLeafRadialGeometry(phi, l, N)

# Pasley (non convex - TODO)
# leaf_rp.la, leaf_rp.lb, leaf_rp.lmax, leaf_rp.ln, leaf_rp.r, leaf_rp.dx = 5.4, 4.3, 5.4 + 4.3 + 0.8, 0.8, 1, 0.5
# phi = np.array([-90, -45, -25, -15, 0, 10, 20, 30, 40, 45, 55, 60, 65, 90]) / 180. * np.pi
# l = np.array([0.8, 0.5, 1.5, 2, 2.6, 2, 2.6, 2.7, 3.1, 3.3, 2.6, 3.3, 3.8, 5.4])
# assert(l.shape == phi.shape)
# N = 200  # N is rather high for testing
# leaf_rp.createLeafRadialGeometry(phi, l, N)

# leaf_rp.la, leaf_rp.lb, leaf_rp.lmax, leaf_rp.ln, leaf_rp.r, leaf_rp.dx = 5, 1, 11, 5, 1, 0.25
# phi = np.array([-90.0, -67.5, -45, -22.5, 0, 22.5, 45, 67.5, 90]) / 180.0 * np.pi
# l = np.array([5.0, 2, 5, 2, 5, 2, 5, 2, 5])
# assert l.shape == phi.shape
# N = 200  # N is rather high for testing
# leaf_rp.createLeafRadialGeometry(phi, l, N)

leaf_rp.areaMax = 50  # cm2

plant.setOrganRandomParameter(leaf_rp)

plant.initialize_empty()

leaf = pb.Leaf(plant, 1, 0, plant.getSeed(), 0)  # create a leaf organ with the given random parameter
leaf.addNode(plant.getSeed().getNode(0), 0, 0.0)  # base node required before plant.simulate(); reuse seed's node ID so SegmentAnalyser indices stay contiguous
plant.getSeed().addChild(leaf)  # defined in Organism

print("**********")
plant.simulate(20)

for o in plant.getOrgans():
    print(o)

vp.plot_plant(plant, "creationTime")

print("fin")
