"""Carrot"""

import numpy as np

import plantbox as pb
import plantbox.visualisation.vtk_plot as vp
from plantbox.visualisation.vtk_animate import AnimateRoots


def carrot_radius(age, z, r0=0.05, rmax=2.5, a_s=30.0, k=6.4e-4, p=1.5, Ls=15, q=2):
    """Secondary-growth radius model r(a, z) = r0 + (rmax - r0) * g(a) * h(z).

    g(a) = 1 - exp(-k * (a - a_s)^p) ramps up secondary growth once age a exceeds a_s.
    h(z) = exp(-(z / Ls)^q) concentrates thickening near the root crown (z=0).

    age: root age [day]
    z:   distance from the root crown [cm]
    r0:  primary-growth (unthickened) radius [cm]
    rmax: fully thickened radius [cm]
    a_s: age at which secondary growth starts [day]
    """
    if age <= a_s:
        return r0
    g = 1.0 - np.exp(-k * (age - a_s) ** p)
    h = np.exp(-((z / Ls) ** q))
    return r0 + (rmax - r0) * g * h


class CarrotRoot(pb.Root):
    """Carrot root subclass that implements a secondary-growth radius model (just for visualization)"""

    def getParameter(self, name, addParams):
        """in case of radius, model parameters need to be passed within addParams"""

        if name in ("radius", "a"):
            if self.param().subType == 1:  # only for taproots
                index = addParams.get("index", 0)
                idx = min(int(index), self.getNumberOfNodes() - 1)
                return carrot_radius(self.getAge(), self.getLength(idx), self.param().a, addParams["rmax"], addParams["a_s"], addParams["k"], addParams["p"], addParams.get("Ls"), addParams["q"])

        return super().getParameter(name, addParams)


class Carrot(pb.Plant):
    """Plant subclass that creates CarrotRoot organs instead of the default Root."""

    def __init__(self, seednum=0):
        super().__init__(seednum)
        self._organs = []  # keeps Python-side CarrotRoot instances alive (otherwise the wrapper is garbage collected)

    def createRoot(self, subType, delay, parent, pni):
        r = CarrotRoot(self, subType, delay, parent, pni)
        self._organs.append(r)
        return r


plant = Carrot()

# Open plant and root parameter from a file
path = "../../modelparameter/structural/rootsystem/"
name = "Daucus_carota"
plant.readParameters(path + name + ".xml")

addParams = {"r0": 0.05, "rmax": 3, "a_s": 30.0, "k": 6.4e-4, "p": 1.5, "Ls": 20, "q": 5}

plant.initialize()

sim_time = 120
dt = 1  # simulation step [day], animation shows one frame per step

# plant.simulate(sim_time)

# Animate growth in an interactive vtk window (close window, or press 'e', to continue)
anim = AnimateRoots(plant, add_params=addParams)
anim.root_name = "age"
anim.start(axis="v")
for i in range(int(sim_time / dt)):
    plant.simulate(dt)
    anim.simtime = (i + 1) * dt
    anim.update()
anim.run()

# Export final result (as vtp)
# plant.write("results/example_plant.vtp")


# ana = pb.SegmentAnalyser(plant, addParams)
# ana.write("results/example_plant_segs.vtp")
# vp.plot_plant(ana, "subType")

# organs = plant.getOrgans()
# taproot = organs[0]
# for i in range(taproot.getNumberOfNodes()):
#     print(f"node {i}: z={taproot.getNode(i).z:.2f}, radius={taproot.getParameter('radius', {'index': i}):.4f}")
