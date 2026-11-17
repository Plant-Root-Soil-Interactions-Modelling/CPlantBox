"""Playing around with Maize parameter"""

import numpy as np
from maize_params import get_internode_lenght, get_leaf_blade_lenght, get_leaf_sheath_length, get_leaf_width

import plantbox as pb
import plantbox.visualisation.vtk_plot as vp
from plantbox.visualisation.vtk_animate import AnimatePlant, AnimateRoots

#
# 1. Seed rotation and orientation, ongoing we need an SRP beta and beta_s for intitial rotation
#
# Grass leaf std beta?
# Harmonize toString outputs
#
# MaizeStem - predef all params, adapt srp for all, sanity checks for leaf
# Brace roots
#


class MaizeStem(pb.Stem):
    """Maize stem subclss"""

    def __init__(self, plant, subType, delay, parent, pni):
        super().__init__(plant, subType, delay, parent, pni)

        self.partialIHeading = pb.Vector3d(-0.5, -0.5, 1)

        srp = self.getStemRandomParameter()
        p = self.param()  # specific stem parameters

        lmax = 197  # total stem length [cm]
        r = lmax / 110  # assuming lmax is reached after 110 days

        # print("Stem random parameters")
        # print(srp)
        # print()
        # print("Stem specific parameters")
        # print(p)

        p.theta = 0  # straight up
        p.lb = 0  # mesocotyl [cm]
        p.la = 0
        ln_ = get_internode_lenght()  # fixed to 19 internodes for now
        sum_ln_ = np.sum(ln_)
        scaled_ln = (ln_ / sum_ln_) * lmax  # scale internode lengths to total stem length of lmax
        p.ln = scaled_ln
        p.r = r

        srp.nodalGrowth = 0
        srp.multDelay = 1
        srp.ldelay = 5
        srp.dx = 0.5
        srp.dxMin = 1.0e-3

        print("Minimal resolution fixed to ", self.getStemRandomParameter().dxMin)
        print(self.param().la)
        print(self.param().lb)
        print(self.param().getK())
        print(self.getInitialHeading())

        # # todo sanity checks ... no tillers, etc.
        # plant = self.getOrganism()
        # seed = plant.getSeed()
        # print(seed)

        print("Successor rules:")
        print("srp.successorWhere", srp.successorWhere)
        print("srp.successorOT", srp.successorOT)
        print("srp.successorST", srp.successorST)
        print("srp.successorP", srp.successorP)
        print("srp.successorNo", srp.successorNo)


class MaizeRoot(pb.Root):
    """Maize root subclass - add brace root behavior"""

    def __init__(self, plant, subType, delay, parent, pni):
        super().__init__(plant, subType, delay, parent, pni)
        p = self.param()  # specific stem parameters
        if p.subType == 5:
            print("crown roots")


class Maize(pb.Plant):
    """Plant subclass that creates MaizeStem organs instead of the default Stem."""

    def __init__(self, seednum=0):
        super().__init__(seednum)
        self.number_of_leafs = 0
        self.ll_ = get_leaf_blade_lenght()
        self.sl_ = get_leaf_sheath_length()
        self.lw_ = get_leaf_width()
        self._organs = []

    def initializeReader(self):
        print("Maize Python initializeReader is called", flush=True)
        super().initializeReader()
        glp = pb.GrassLeafRandomParameter(self)
        glp.subType = 0
        self.setOrganRandomParameter(glp)

    def createLeaf(self, subType, delay, parent, pni):
        """grass leaf instead of the default leaf, hard coded parameter values"""
        leaf = pb.GrassLeaf(self, subType, delay, parent, pni)
        lrp = leaf.getGrassLeafRandomParameter()
        lp = leaf.param()
        # hack to set parameters
        # print("Leaf specific parameters")
        # print(leaf.param())
        # a	0.0297928
        # bladeAngle	0.3
        # sheathLength	5
        # leafGrowthDuration	20
        lp.bladeLength = self.ll_[self.number_of_leafs]  # * 0.0
        lp.sheathLength = self.sl_[self.number_of_leafs] * 0.0
        lp.bladeWidth = self.lw_[self.number_of_leafs]
        lp.bladeBending = 0.01
        self.number_of_leafs += 1
        self._organs.append(leaf)
        return leaf

    def createStem(self, subType, delay, parent, pni):
        """simplified stem with hard coded internode distances"""
        r = MaizeStem(self, subType, delay, parent, pni)
        self._organs.append(r)
        return r

    def createRoot(self, subType, delay, parent, pni):
        """prototyping brace root behaviour"""
        r = MaizeRoot(self, subType, delay, parent, pni)
        self._organs.append(r)
        return r


plant = Maize()

# Open plant and root parameter from a file
path = "../../modelparameter/structural/plant/"
name = "2020-maize"
plant.readParameters(path + name + ".xml")

#
# Seed
#
srp = plant.getOrganRandomParameter(pb.seed)
srp[0].seedPos.x = 0
srp[0].seedPos.y = 0

# Delay timing type
srp[0].delayDefinitionRoot = pb.dd_distance
srp[0].delayDefinitionShoot = pb.dd_time_lat

# Shoot borne definition
srp[0].nC = 7
srp[0].firstSB = 0
srp[0].delaySB = 3

#
# Main stem
#
# we want to add (for now)
# <parameter name="successor" ruleId="1" numLat="nC" subType="5" organType ="2"  where="0" percentage="1"/>
strp = plant.getOrganRandomParameter(pb.stem)

print("srp.successorWhere", strp[1].successorWhere)
print("srp.successorOT", strp[1].successorOT)
print("srp.successorST", strp[1].successorST)
print("srp.successorP", strp[1].successorP)
print("srp.successorNo", strp[1].successorNo)

print("***")
successorWhere = strp[1].successorWhere
successorOT = strp[1].successorOT
successorST = strp[1].successorST
successorP = strp[1].successorP
successorNo = strp[1].successorNo

successorWhere.append([0])
successorOT.append([pb.root])
successorST.append([5])
successorP.append([1])
successorNo.append(int(srp[0].nC))

strp[1].successorWhere = successorWhere
strp[1].successorOT = successorOT
strp[1].successorST = successorST
strp[1].successorP = successorP
strp[1].successorNo = successorNo


print("srp.successorWhere", strp[1].successorWhere)
print("srp.successorOT", strp[1].successorOT)
print("srp.successorST", strp[1].successorST)
print("srp.successorP", strp[1].successorP)
print("srp.successorNo", strp[1].successorNo)

#
# Crown root
#
rrp = plant.getOrganRandomParameter(pb.root)
rrp[5].ldelays = srp[0].firstSB
rrp[5].ldelays = 0.0
rrp[5].multDelay = srp[0].delaySB


#  double ldelays = 0.;        ///< Standard deviation of the lateral emergence delay [day]
#     double multDelay = 0;

plant.initialize()
plant.writeParameters(name + "_modified" + ".xml")


sim_time = 120
dt = 1  # simulation step [day], animation shows one frame per step
N = sim_time / dt

# for i in range(int(N)):
#     plant.simulate(dt)
# ana = pb.SegmentAnalyser(plant)
# ana.addAge(sim_time)
# vp.plot_roots(ana, "age")  # Interactive plot, using vtk


# Animate growth in an interactive vtk window (close window, or press 'e', to continue)
anim = AnimatePlant(plant)
anim.param_name = "subType"  # "creationTime"
anim.start(axis="v")
anim.plant = True
anim.avi_name = "maize_"
for i in range(int(sim_time / dt)):
    print(i)
    plant.simulate(dt)
    anim.simtime = (i + 1) * dt
    anim.update()
anim.run()
anim.make_video(output_file="maize.mp4", fps=20, codec="libx264")
