"""
Minimal GrassLeaf simulation from XML
======================================
Loads grassleaf_parameters.xml, simulates 25 days, plots via plot_plant.

Run from this directory:
    python grassleaf_xml.py
"""

import os

import numpy as np

import plantbox as pb
import plantbox.visualisation.vtk_plot as vp


class Poaceae(pb.Plant):
    """Plant subclass that registers GrassLeafRandomParameter as the leaf prototype."""

    def initializeReader(self):
        super().initializeReader()
        gl_proto = pb.GrassLeafRandomParameter(self)
        gl_proto.subType = 0
        self.setOrganRandomParameter(gl_proto)

    def createLeaf(self, subType, delay, parent, pni):
        return pb.GrassLeaf(self, subType, delay, parent, pni)


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
