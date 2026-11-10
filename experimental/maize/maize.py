"""Carrot"""

import numpy as np

import plantbox as pb
import plantbox.visualisation.vtk_plot as vp
from plantbox.visualisation.vtk_animate import AnimateRoots


class MaizeStem(pb.Stem):
    """Maize stem subclss"""

    def __init__(self, plant, subType, delay, parent, pni):
        super().__init__(plant, subType, delay, parent, pni)

        print("Stem random parameters")
        print(self.getStemRandomParameter())
        print()
        print("Stem specific parameters")
        print(self.param())


class Maize(pb.Plant):
    """Plant subclass that creates MaizeStem organs instead of the default Stem."""

    def __init__(self, seednum=0):
        super().__init__(seednum)
        self._organs = []  # keeps Python-side CarrotRoot instances alive (otherwise the wrapper is garbage collected)

    def createStem(self, subType, delay, parent, pni):
        r = MaizeStem(self, subType, delay, parent, pni)
        self._organs.append(r)
        return r


plant = Maize()

# Open plant and root parameter from a file
path = "../../modelparameter/structural/plant/"
name = "2020-maize"
plant.readParameters(path + name + ".xml")

srpm = plant.getOrganRandomParameter(pb.seed)
srpm[0].seedPos.x = 0
srpm[0].seedPos.y = 0
# print(srpm[0])

plant.initialize()

sim_time = 20
dt = 1  # simulation step [day], animation shows one frame per step

plant.simulate(sim_time)

# # Animate growth in an interactive vtk window (close window, or press 'e', to continue)
# anim = AnimateRoots(plant)
# anim.root_name = "age"
# anim.start(axis="v")
# for i in range(int(sim_time / dt)):
#     plant.simulate(dt)
#     anim.simtime = (i + 1) * dt
#     anim.update()
# anim.run()


# Interactive plot, using vtk
vp.plot_plant(plant, "subType")


# organs = plant.getOrgans()
# taproot = organs[0]
# for i in range(taproot.getNumberOfNodes()):
#     print(f"node {i}: z={taproot.getNode(i).z:.2f}, radius={taproot.getParameter('radius', {'index': i}):.4f}")
