"""Playing around with Maize parameter"""

import numpy as np

import plantbox as pb
import plantbox.visualisation.vtk_plot as vp
from plantbox.visualisation.vtk_animate import AnimateRoots


def get_internode_lenght():
    internode_length = np.array(
        [
            0.02,
            0.02,
            0.02,
            0.02,
            0.02,  # 1–5: basal
            0.085,
            0.085,
            0.085,
            0.085,  # 6–9: lower-middle
            0.16,
            0.16,
            0.16,
            0.16,  # 10–13: middle
            0.115,
            0.115,
            0.115,
            0.115,  # 14–17: upper
            0.02,
            0.02,  # 18–19: below tassel
        ]
    )  # [m]
    return internode_length * 100  # [m] -> [cm]


class MaizeStem(pb.Stem):
    """Maize stem subclss"""

    def __init__(self, plant, subType, delay, parent, pni):
        super().__init__(plant, subType, delay, parent, pni)

        lmax = 197  # total stem length [cm]
        r = lmax / 110  # assuming lmax is reached after 110 days

        print("Stem random parameters")
        print(self.getStemRandomParameter())
        print()
        print("Stem specific parameters")

        print(self.param())
        self.param().theta = 0  # straight up
        self.param().lb = 2  # mesocotyl [cm]
        self.param().la = 0
        ln_ = get_internode_lenght()  # fixed to 19 internodes for now
        sum_ln_ = np.sum(ln_)
        scaled_ln = (ln_ / sum_ln_) * lmax  # scale internode lengths to total stem length of lmax
        self.param().ln = scaled_ln
        self.param().r = 0.1  # not working, why? callbacks?

        # todo sanity checks ... no tillers, etc.


class Maize(pb.Plant):
    """Plant subclass that creates MaizeStem organs instead of the default Stem."""

    def __init__(self, seednum=0):
        super().__init__(seednum)
        self._organs = []  # keeps Python-side CarrotRoot instances alive (otherwise the wrapper is garbage collected)

    def createLeaf(self, subType, delay, parent, pni):
        leaf = pb.GrassLeaf(self, subType, delay, parent, pni)
        # hack to set area
        print("Leaf specific parameters")
        print(leaf.param())
        self._organs.append(leaf)
        return leaf

    def initializeReader(self):
        print("Maize Python initializeReader is called", flush=True)
        super().initializeReader()
        glp = pb.GrassLeafRandomParameter(self)
        glp.subType = 0
        self.setOrganRandomParameter(glp)

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

plant.writeParameters(name + "_modified" + ".xml")


# Interactive plot, using vtk
vp.plot_plant(plant, "subType")
