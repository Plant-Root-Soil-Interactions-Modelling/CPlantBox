"""How to define a custom Plant subclass for Poaceae that creates GrassLeaf organs instead of the default Leaf"""

import math

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

import plantbox as pb
import plantbox.visualisation.vtk_plot as vp


class Poaceae(pb.Plant):
    """Plant subclass that creates GrassLeaf organs instead of the default Leaf."""

    def createLeaf(self, subType, delay, parent, pni):
        return pb.GrassLeaf(self, subType, delay, parent, pni)

    def initializeReader(self):
        super().initializeReader()
        glp = pb.GrassLeafRandomParameter(self)
        glp.subType = 0
        self.setOrganRandomParameter(glp)


# Simulate a Poaceae plant with GrassLeaf organs
plant = Poaceae()

xml_path = os.path.join(os.path.dirname(__file__), "grassleaf_parameters.xml")
plant.readParameters(xml_path)
plant.initialize(verbose=False)

total_days = 25.0
dt = 0.5
steps = int(total_days / dt)

for i in range(steps):
    plant.simulate(dt, verbose=False)

ana = pb.SegmentAnalyser(plant)
ana.addAge(total_days)
vp.plot_plant(ana, "age")
