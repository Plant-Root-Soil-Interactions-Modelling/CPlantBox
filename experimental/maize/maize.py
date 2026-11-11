"""Playing around with Maize parameter"""

import numpy as np
from maize_params import get_internode_lenght, get_leaf_blade_lenght, get_leaf_sheath_length, get_leaf_width

import plantbox as pb
import plantbox.visualisation.vtk_plot as vp
from plantbox.visualisation.vtk_animate import AnimatePlant, AnimateRoots


class MaizeStem(pb.Stem):
    """Maize stem subclss"""

    def __init__(self, plant, subType, delay, parent, pni):
        super().__init__(plant, subType, delay, parent, pni)

        srp = self.getStemRandomParameter()
        p = self.param()  # specific stem parameters

        lmax = 197  # total stem length [cm]
        r = lmax / 110  # assuming lmax is reached after 110 days

        print("Stem random parameters")
        print(srp)
        print()
        print("Stem specific parameters")
        print(p)

        p.theta = 0  # straight up
        p.lb = 2  # mesocotyl [cm]
        p.la = 0
        ln_ = get_internode_lenght()  # fixed to 19 internodes for now
        sum_ln_ = np.sum(ln_)
        scaled_ln = (ln_ / sum_ln_) * lmax  # scale internode lengths to total stem length of lmax
        p.ln = scaled_ln
        p.r = r

        srp.multDelay = 1
        srp.ldelay = 5

        # # todo sanity checks ... no tillers, etc.
        # plant = self.getOrganism()
        # seed = plant.getSeed()
        # print(seed)


class Maize(pb.Plant):
    """Plant subclass that creates MaizeStem organs instead of the default Stem."""

    def __init__(self, seednum=0):
        super().__init__(seednum)
        self.number_of_leafs = 0
        self.ll_ = get_leaf_blade_lenght()
        self.sl_ = get_leaf_sheath_length()
        self.lw_ = get_leaf_width()
        self._organs = []  # keeps Python-side CarrotRoot instances alive (otherwise the wrapper is garbage collected)

    def createLeaf(self, subType, delay, parent, pni):
        leaf = pb.GrassLeaf(self, subType, delay, parent, pni)
        lp = leaf.param()
        # hack to set parameters
        # print("Leaf specific parameters")
        # print(leaf.param())
        # a	0.0297928
        # bladeAngle	0.3
        # sheathLength	5
        #
        # leafGrowthDuration	20
        lp.bladeLength = self.ll_[self.number_of_leafs]
        lp.sheathLength = self.sl_[self.number_of_leafs]
        lp.bladeWidth = self.lw_[self.number_of_leafs]
        lp.bladeBending = 0.01
        self.number_of_leafs += 1
        print(self.number_of_leafs)
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
print(srpm[0])
srpm[0].delayDefinitionRoot = pb.dd_distance
srpm[0].delayDefinitionShoot = pb.dd_time_lat


plant.initialize()

sim_time = 120
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
