import numpy as np
from maize_params import *

import plantbox as pb
import plantbox.visualisation.vtk_plot as vp
from plantbox.visualisation.vtk_animate import AnimatePlant, AnimateRoots


class MaizeStem(pb.Stem):
    """Maize stem"""

    def __init__(self, plant, subType, delay, parent, pni):
        super().__init__(plant, subType, delay, parent, pni)

        self.number_of_roots = {}  # at each stem node index (needed by MaizeRoot)

        srp = self.getStemRandomParameter()
        p = self.param()  # specific stem parameters

        if p.subType == 1:  # only for mainstem

            lmax = 197  # total stem length [cm]
            r = lmax / 110  # assuming lmax is reached after 110 days

            # print("Stem random parameters")
            # print(srp)
            # print()
            print("Stem specific parameters")
            print(p)
            p.theta = 0  # straight up
            p.lb = 2.0  # mesocotyl [cm]
            p.la = 0.0
            ln_ = get_internode_lenght()  # fixed to 19 internodes for now
            sum_ln_ = np.sum(ln_)
            scaled_ln = (ln_ / sum_ln_) * lmax  # scale internode lengths to total stem length of lmax
            p.ln = scaled_ln
            p.r = r
            p.ldelay = 0  # unused, ldelay are set for the leafs and shoot borne roots (since srp[0].delayDefinitionShoot = pb.dd_time_self)

            srp.nodalGrowth = 0  # 0 = sequentiall growth, 1=simultaneous growth
            srp.multDelay = 0  # unused (since srp[0].delayDefinitionShoot = pb.dd_time_self)
            srp.dx = 0.5
            srp.dxMin = 1.0e-3

            print("Minimal resolution fixed to ", self.getStemRandomParameter().dxMin)
            print(self.param().la)
            print(self.param().lb)
            print(self.param().getK())
            print(self.getInitialHeading())

            print("Successor rules:")
            print("srp.successorWhere", srp.successorWhere)
            print("srp.successorOT", srp.successorOT)
            print("srp.successorST", srp.successorST)
            print("srp.successorP", srp.successorP)
            print("srp.successorNo", srp.successorNo)
            print("*")

            # # todo sanity checks ... no tillers, etc.
            # plant = self.getOrganism()
            # seed = plant.getSeed()
            # print(seed)

    def simulate(self, dt, verbose=False):
        super().simulate(dt, verbose)
        # Make stem radius age dependent
        p = self.param()  # specific stem parameters
        days_after_sowing = self.getPlant().getSimTime()
        min_a = 0.2
        max_a = 1.5
        t = days_after_sowing / 110  # assuming lmax is reached after 110 days
        t = min(1.0, t)
        p.a = min_a * (1 - t) + max_a * t


class MaizeRoot(pb.Root):
    """Maize root subclass - add brace root behavior"""

    def __init__(self, plant, subType, delay, parent, pni):
        super().__init__(plant, subType, delay, parent, pni)
        p = self.param()  # specific stem parameters
        if p.subType == 5:
            rrp = self.getRootRandomParameter()
            parent = self.getParent()
            if pni not in parent.number_of_roots:
                parent.number_of_roots[pni] = 0
            parent.number_of_roots[pni] += 1
            rrp.multDelay = 0.0  # dont' use it
            lni = parent.getCreatedLinkingNode()  # current linking node index = whorl number
            print("Crown root", parent.number_of_roots[pni], "linking node number", lni)
            nC = get_crown_root_number()[lni]
            n = parent.number_of_roots[pni]
            whorl_beta = 2 * np.pi / 6.0 * lni  # no idea how to choose it, but propably should not be equal to all whorls
            beta = whorl_beta + 2 * np.pi * n / nC
            theta = np.pi / 2.0 - get_crown_root_angles()[lni] / 180.0 * np.pi
            self.setPartialInitialHeading(theta, beta)
            # todo: make sure there are is no soil geometry attached to the tropisms


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
        lni = parent.getCreatedLinkingNode()  # current linking node index = phytomere number
        delay = get_leaf_emergence_days()[lni]
        leaf = pb.GrassLeaf(self, subType, delay, parent, pni)
        lrp = leaf.getGrassLeafRandomParameter()
        lrp.rotDev = 0.2  # add randomness to the initial leaf rotation
        lp = leaf.param()
        lp.bladeLength = self.ll_[self.number_of_leafs]  # * 0.0
        lp.sheathLength = self.sl_[self.number_of_leafs]  # * 0.0
        lp.bladeWidth = self.lw_[self.number_of_leafs]
        lp.bladeBending = 0.01
        # todo
        self.number_of_leafs += 1
        self._organs.append(leaf)
        return leaf

    def createStem(self, subType, delay, parent, pni):
        """simplified stem with hard coded internode distances"""
        print("Creating MaizeStem with subType:", subType, "and delay:", delay)
        delay = 2  # days
        r = MaizeStem(self, subType, delay, parent, pni)
        self._organs.append(r)
        return r

    def createRoot(self, subType, delay, parent, pni):
        """prototyping brace root behaviour"""
        if subType == 5:  # in case of brace root fix delay time based on data per whorl number
            lni = parent.getCreatedLinkingNode()  # current linking node index = whorl number
            cre = get_crown_root_emergence()
            delay = cre[lni]
        r = MaizeRoot(self, subType, delay, parent, pni)
        self._organs.append(r)
        return r


plant = Maize()

# Open plant and root parameter from a file
# path = "../../modelparameter/structural/plant/"
# name = "2020-maize"
path = ""
# name = "2020-maize_modified"
name = "Zeamays_synMRI"
plant.readParameters(path + name + ".xml")

#
# Seed
#

srp = plant.getOrganRandomParameter(pb.seed)
srp[0].seedPos.x = 0
srp[0].seedPos.y = 0
srp[0].seedPos.z = -3

# Delay timing type
srp[0].delayDefinitionRoot = pb.dd_distance
srp[0].delayDefinitionShoot = pb.dd_time_self  # each lateral defines it's own delay

# Shoot borne definition
srp[0].nC = 0  # turn it off
srp[0].delaySB = 7
srp[0].firstSB = 1e9  # turn it off

#
# Main stem
#
# we want to add (for now)
# <parameter name="successor" ruleId="1" numLat="nC" subType="5" organType ="2"  where="0" percentage="1"/>

strp = plant.getOrganRandomParameter(pb.stem)

print(strp[1].gf)  # Growth function (1=negative exponential, 2=linear)


# successorWhere = strp[1].successorWhere
# successorOT = strp[1].successorOT
# successorST = strp[1].successorST
# successorP = strp[1].successorP
# successorNo = strp[1].successorNo

# nC = get_crown_root_number()
# for i in range(0, 6):
#     successorWhere.append([int(i)])  # whorl index; check at which node the coleoptilar node is; that's the starting node
#     successorOT.append([pb.root])
#     successorST.append([5])
#     successorP.append([1])
#     successorNo.append(nC[i])  # inceasing per root whorle (approx int(nC + 2 * i) for nC = 5)

# strp[1].successorWhere = successorWhere
# strp[1].successorOT = successorOT
# strp[1].successorST = successorST
# strp[1].successorP = successorP
# strp[1].successorNo = successorNo

# print("srp.successorWhere", strp[1].successorWhere)
# print("srp.successorOT", strp[1].successorOT)
# print("srp.successorST", strp[1].successorST)
# print("srp.successorP", strp[1].successorP)
# print("srp.successorNo", strp[1].successorNo)


#
# Crown root
#
rrp = plant.getOrganRandomParameter(pb.root)
rrp[5].ldelay = 0.0
rrp[5].ldelays = 0.0
rrp[5].multDelay = 0  # roots emerge simultaneously foreach root whorl
switch_age = 30  # for example, switch from exotropism to gravitropism after 5 days
exogravitropism = ExoGraviTropism(plant, rrp[5].tropismN, rrp[5].tropismS, switch_age)

#
# Leaf
#
lrp = plant.getOrganRandomParameter(pb.leaf)
print("Leaf", lrp[1].subType)
print("Leaf delay", lrp[1].ldelay)
print("Leaf multDelay", lrp[1].multDelay)


plant.initialize()  # initialize sets the tropisms
plant.setTropism(exogravitropism, pb.root, 5)
plant.writeParameters(name + "_modified" + ".xml")


sim_time = 120  # 120 # days
anim_time = 12  # seconds
fps = 20  # frames per second
dt = sim_time / (anim_time * fps)  # simulation step [day], animation shows one frame per step
N = int(sim_time / dt)


# for i in range(int(N)):
#     plant.simulate(dt)
# ana = pb.SegmentAnalyser(plant)
# ana.addAge(sim_time)
# vp.plot_plant(ana, "subType")  # Interactive plot, using vtk "age"


insert_times = [8, 12, 25, 35, 45, 50, 60, 65, 70, 75, 80]
text = ["VE", "V2", "V4", "V6", "V8", "V10", "V12", "V14", "V16–V18", "VT", "R1"]
text_dec = [
    "Internodes remain very short",
    "Initial stem axis development",
    "Early internode differentiation",
    "Start of rapid basal internode elongation",
    "Basal internodes elongate rapidly",
    "Strong elongation of basal and middle internodes",
    "Middle internodes elongate rapidly",
    "Upper internodes develop; basal internodes approach final length",
    "Upper internodes elongate",
    "Plant approaches final height; most internodes fully elongated",
    "Internode elongation essentially complete",
]

# Animate growth in an interactive vtk window (close window, or press 'e', to continue)
anim = AnimatePlant(plant)
anim.param_name = "organType"  # "creationTime"
anim.set_captions(insert_times, text, text_dec)
anim.show_time = True
anim.start(axis="v")
anim.plant = True
anim.avi_name = "maize_"
for i in range(int(sim_time / dt)):
    plant.simulate(dt)
    anim.simtime = (i + 1) * dt
    anim.update()
anim.run()
anim.make_video(output_file="maize.mp4", fps=fps, codec="libx264")
