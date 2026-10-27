import sys; sys.path.append(".."); sys.path.append("../src/")
import unittest

import plantbox as pb

import numpy as np

path = "../modelparameter/structural/plant/"


class TestMycorrhizalPlant(unittest.TestCase):

    def myc_plant(self, lmbd = 2., hed = 1., distTH = 0., gravitropic = False, seed = 1):
        """ a mycorrhizal plant with a single colonized tap root, hyphae emerge from the colonized root nodes and branch """
        self.plant = pb.MycorrhizalPlant(seed)
        srp = pb.SeedRandomParameter(self.plant)
        srp.subType = 0
        self.plant.setOrganRandomParameter(srp)
        self.rrp = pb.MycorrhizalRootRandomParameter(self.plant)
        self.rrp.name, self.rrp.subType, self.rrp.lmax, self.rrp.r, self.rrp.dx = "taproot", 1, 10., 2., 0.5
        self.rrp.lmbd, self.rrp.highresolution, self.rrp.hyphalEmergenceDensity = lmbd, 0., hed
        if gravitropic:  # straight down
            self.rrp.tropismT, self.rrp.tropismS = 1, 0.
        self.plant.setOrganRandomParameter(self.rrp)
        self.hrp = pb.HyphaeRandomParameter(self.plant)  # root hyphae are always of sub type 1
        self.hrp.name, self.hrp.subType, self.hrp.dx = "hyphae", 1, 0.05
        self.hrp.v, self.hrp.vs, self.hrp.lmax, self.hrp.lmaxs, self.hrp.la, self.hrp.lb, self.hrp.ln = 1., 0., 1., 0., 0.2, 0.2, 0.2
        self.hrp.distTH = distTH
        self.hrp.successor = [[1]]
        self.hrp.successorP = [[1.]]
        self.plant.setOrganRandomParameter(self.hrp)
        self.plant.initialize(False)

    def xml_plant(self, name = "Heliantus_Pagès_2013.xml"):
        """ a plant from a parameter file, hyphae parameters are added (the file has none) """
        self.plant = pb.MycorrhizalPlant(1)
        self.plant.readParameters(path + name, fromFile = True, verbose = False)
        hrp = pb.HyphaeRandomParameter(self.plant)
        hrp.subType = 1
        self.plant.setOrganRandomParameter(hrp)
        self.plant.initialize(False)

    def simulate(self, time, dt = 0.25):
        for i in range(0, int(round(time / dt))):
            self.plant.simulate(dt, False)

    def test_read_parameters(self):
        """ tests that a parameter file is read into mycorrhizal parameter sets """
        p = pb.MycorrhizalPlant()
        p.readParameters(path + "Heliantus_Pagès_2013.xml", fromFile = True, verbose = False)
        seeds = p.getOrganRandomParameter(pb.OrganTypes.seed)
        roots = p.getOrganRandomParameter(pb.OrganTypes.root)
        stems = p.getOrganRandomParameter(pb.OrganTypes.stem)
        leafs = p.getOrganRandomParameter(pb.OrganTypes.leaf)
        self.assertEqual([len(seeds), len(roots[1:]), len(stems[1:]), len(leafs[1:])], [1, 3, 3, 1],
                         "read parameters: read wrong number of random parameter from xml")
        for rp in roots:
            self.assertEqual(rp.__class__.__name__, "MycorrhizalRootRandomParameter", "read parameters: root parameters are not mycorrhizal")
        self.assertEqual(len(p.getOrganRandomParameter(pb.hyphae)), 1, "read parameters: expected the default hyphae parameter set")

    def test_simulate_xml(self):
        """ tests the simulation of a plant from a parameter file, all roots are mycorrhizal roots """
        self.xml_plant()
        self.simulate(10., 1.)
        roots = self.plant.getOrgans(pb.root)
        self.assertGreater(len(roots), 1, "simulate: no laterals")
        for r in roots:
            self.assertEqual(r.__class__.__name__, "MycorrhizalRoot", "createRoot: root is not a MycorrhizalRoot")
        self.assertGreater(len(self.plant.getOrgans(pb.hyphae)), 0, "simulate: no hyphae emerged")
        self.assertGreater(len(self.plant.getOrgans(pb.stem)), 0, "simulate: no stem")
        ana = pb.SegmentAnalyser(self.plant)
        ana.addData("colonization", self.plant.getNodeColonizations(pb.root))
        ana.addData("anastomosis", self.plant.getAnastomosisPoints(pb.hyphae))
        ana.write("mycorrhizalplant.vtp", ["radius", "subType", "creationTime", "organType", "colonization", "anastomosis"])

    def test_missing_laterals(self):
        """ tests that probabilistic branching without laterals at a branching point does not crash """
        self.xml_plant("test_missing_laterals.xml")
        self.simulate(20., 1.)

    def node_colonizations(self):
        self.xml_plant()
        self.simulate(10., 1.)
        return self.plant.getOrgans(pb.root), np.array(self.plant.getNodeColonizations(pb.root)), np.array(self.plant.getNodeColonizationTime(pb.root))

    def test_node_colonizations(self):
        """ tests getNodeColonizations() and getNodeColonizationTime() against the colonization of the roots (without the root tips) """
        roots, col, ct = self.node_colonizations()
        self.assertEqual(len(col), self.plant.getNumberOfNodes(), "getNodeColonizations: unexpected length")
        self.assertEqual(len(ct), self.plant.getNumberOfNodes(), "getNodeColonizationTime: unexpected length")
        self.assertTrue(set(np.unique(col)) <= {0, 1, 2, 3}, "getNodeColonizations: unexpected colonization state")
        self.assertGreater(np.sum(col > 0), 0, "getNodeColonizations: nothing colonized")
        for r in roots:
            for i in range(1, r.getNumberOfNodes() - 1):  # node 0 is the parent node (or the seed node), it is set by the parent
                nid = r.getNodeId(i)
                self.assertEqual(col[nid], r.getNodeColonization(i), "getNodeColonizations: unexpected state")
                self.assertEqual(ct[nid], r.getNodeColonizationTime(i), "getNodeColonizationTime: unexpected time")

    def test_node_colonizations_tips(self):
        """ tests getNodeColonizations() at the root tips """
        roots, col, ct = self.node_colonizations()
        tips = [(r.getNodeId(r.getNumberOfNodes() - 1), r.getNodeColonization(r.getNumberOfNodes() - 1)) for r in roots]
        self.assertGreater(sum(1 for _, c in tips if c > 0), 0, "getNodeColonizations: no colonized root tip in this example")
        for nid, c in tips:
            self.assertEqual(col[nid], c, "getNodeColonizations: unexpected state at the root tip")

    def test_anastomosis_points(self):
        """ tests getAnastomosisPoints(): 1 at the last node of each merged hypha """
        self.myc_plant(hed = 3., distTH = 0.3)
        self.simulate(3.)
        hyphae = self.plant.getOrgans(pb.hyphae)
        ap = np.array(self.plant.getAnastomosisPoints(pb.hyphae))
        self.assertEqual(len(ap), self.plant.getNumberOfNodes(), "getAnastomosisPoints: unexpected length")
        merged = {h.getNodeId(h.getNumberOfNodes() - 1) for h in hyphae if h.getParameter("mergePointID") > -1}
        self.assertGreater(len(merged), 0, "getAnastomosisPoints: no anastomosis in this example")
        self.assertEqual(set(np.where(ap == 1)[0]), merged, "getAnastomosisPoints: unexpected anastomosis points")

    def test_node_tips(self):
        """ tests getNodeTips(): 1 at the last node of each organ of the given type (and of the base organs, i.e. the seed) """
        self.myc_plant()
        self.simulate(3.)
        for ot in [pb.root, pb.hyphae]:
            tips = np.array(self.plant.getNodeTips(ot))
            self.assertEqual(len(tips), self.plant.getNumberOfNodes(), "getNodeTips: unexpected length")
            expected = {o.getNodeId(o.getNumberOfNodes() - 1) for o in self.plant.getOrgans(ot)}
            seed = self.plant.getSeed()
            expected.add(seed.getNodeId(seed.getNumberOfNodes() - 1))
            self.assertEqual(set(np.where(tips == 1)[0]), expected, "getNodeTips: unexpected tips")

    def test_hyphal_tree_index(self):
        """ tests that the plant numbers the hyphal trees of the root hyphae consecutively """
        self.myc_plant()
        self.simulate(3.)
        indices = sorted(int(h.getParameter("hyphalTreeIndex")) for h in self.plant.getOrgans(pb.hyphae) if h.getParent().organType() == pb.root)
        self.assertGreater(len(indices), 0, "hyphal tree index: no root hyphae")
        self.assertEqual(indices, list(range(len(indices))), "hyphal tree index: root hyphae are not numbered consecutively")

    def test_simulate_hyphae(self):
        """ tests that simulateHyphae() only grows the hyphae """
        self.myc_plant()
        self.simulate(2.)
        root_length = [r.getLength() for r in self.plant.getOrgans(pb.root)]
        hyphae_length = sum(h.getLength() for h in self.plant.getOrgans(pb.hyphae))
        simtime = self.plant.getSimTime()
        self.plant.simulateHyphae(0.25, False)
        self.assertEqual([r.getLength() for r in self.plant.getOrgans(pb.root)], root_length, "simulateHyphae: roots grew")
        self.assertGreater(sum(h.getLength() for h in self.plant.getOrgans(pb.hyphae)), hyphae_length, "simulateHyphae: hyphae did not grow")
        self.assertAlmostEqual(self.plant.getSimTime(), simtime + 0.25, 12, "simulateHyphae: simulation time not increased")

    def test_simulate_colonization(self):
        """ tests that simulateColonization() colonizes the roots and simulateHyphalGrowth() creates hyphae, without root growth """
        self.myc_plant(lmbd = 0.5, hed = 4.)  # every colonized segment (dx = 0.5) gets hyphae
        self.simulate(1.)
        roots = self.plant.getOrgans(pb.root)
        root_length = [r.getLength() for r in roots]

        def root_hyphae():  # new hyphae have age 0 and a single node, getOrgans() skips them
            return sum(1 for r in roots for i in range(r.getNumberOfChildren()) if r.getChild(i).organType() == pb.hyphae)

        n_col = np.sum(np.array(self.plant.getNodeColonizations(pb.root)) > 0)
        n_hyphae = root_hyphae()
        for i in range(10):
            self.plant.simulateColonization(0.25, False)
        self.assertEqual([r.getLength() for r in self.plant.getOrgans(pb.root)], root_length, "simulateColonization: roots grew")
        self.assertGreater(np.sum(np.array(self.plant.getNodeColonizations(pb.root)) > 0), n_col, "simulateColonization: no new colonization")
        self.assertEqual(root_hyphae(), n_hyphae, "simulateColonization: hyphae emerged")
        self.plant.simulateHyphalGrowth(0.25, False)
        self.assertGreater(root_hyphae(), n_hyphae, "simulateHyphalGrowth: no new hyphae")

    def test_change_geometry(self):
        """ tests that changeGeometry() confines the hyphae to the given geometry """
        box = pb.SDF_PlantBox(1., 1., 20.)
        outside = []
        for geometry in [False, True]:
            self.myc_plant(gravitropic = True)
            if geometry:
                self.plant.changeGeometry(pb.hyphae, box)
            self.simulate(3.)
            nodes = [n for h in self.plant.getOrgans(pb.hyphae) for n in h.getNodes()]
            outside.append(sum(1 for n in nodes if box.getDist(n) > 1.e-6))
        self.assertGreater(outside[0], 0, "changeGeometry: example hyphae do not leave the box without geometry")
        self.assertEqual(outside[1], 0, "changeGeometry: hyphae nodes outside the geometry")

    def test_turn_off_side_petri_dish(self):
        """ tests that turnOffSidePetriDish() deactivates organs without a node in the given region, and that they stop growing """
        self.myc_plant()
        self.simulate(2.)
        organs = self.plant.getOrgans()
        x0 = np.median([max(n.x for n in o.getNodes()) for o in organs])  # about half of the organs reach into the region x > x0
        self.plant.turnOffSidePetriDish(x0, 0., 100., -100.)
        in_region = [any(n.x > x0 and n.z < 0. and -100. < n.y < 100. for n in o.getNodes()) for o in organs]
        self.assertGreater(sum(in_region), 0, "turnOffSidePetriDish: no organ in the region")
        self.assertLess(sum(in_region), len(organs), "turnOffSidePetriDish: all organs in the region")
        for o, inside in zip(organs, in_region):
            self.assertEqual(o.isActive(), inside, "turnOffSidePetriDish: unexpected active state")
        hyphae = [o for o, inside in zip(organs, in_region) if o.organType() == pb.hyphae and not inside]
        length = [h.getLength() for h in hyphae]
        self.simulate(0.5)
        self.assertEqual([h.getLength() for h in hyphae], length, "turnOffSidePetriDish: deactivated hyphae grew")

    def test_colonization_rate_laterals(self):
        """ tests that the primary colonization rate of laterals equals that of the tap root (laterals are colonized once per time step) """
        lmbd, dt, time = 0.2, 0.1, 3.
        stats = {1: [0, 0], 2: [0, 0]}  # sub type: [uncolonized nodes, colonized nodes]
        for seed in range(1, 21):
            p = pb.MycorrhizalPlant(seed)
            srp = pb.SeedRandomParameter(p)
            srp.subType = 0
            p.setOrganRandomParameter(srp)
            p0 = pb.MycorrhizalRootRandomParameter(p)  # tap root and laterals reach their final length within the first time steps
            p0.name, p0.subType, p0.lmax, p0.r, p0.dx, p0.gf, p0.la, p0.lb, p0.ln = "taproot", 1, 10., 100., 0.5, 2, 1., 1., 2.
            p0.successor = [[2]]
            p0.successorP = [[1.]]
            p1 = pb.MycorrhizalRootRandomParameter(p)
            p1.name, p1.subType, p1.lmax, p1.r, p1.dx, p1.gf = "lateral", 2, 5., 100., 0.5, 2
            for q in [p0, p1]:
                q.lmbd, q.vi, q.maxAge, q.highresolution, q.hyphalEmergenceDensity = lmbd, 0., 1.e6, 0., 0.
                p.setOrganRandomParameter(q)
            hrp = pb.HyphaeRandomParameter(p)
            hrp.subType = 1
            p.setOrganRandomParameter(hrp)
            p.initialize(False)
            p.simulate(0.2, False)
            before = {r.getId(): [r.getNodeColonization(i) for i in range(r.getNumberOfNodes())] for r in p.getOrgans(pb.root)}
            for i in range(int(round(time / dt))):
                p.simulate(dt, False)
            for r in p.getOrgans(pb.root):
                st = int(r.getParameter("subType"))
                c0 = before.get(r.getId(), [0] * r.getNumberOfNodes())
                for i in range(1, r.getNumberOfNodes()):
                    if i < len(c0) and c0[i] == 0:
                        stats[st][0] += 1
                        stats[st][1] += r.getNodeColonization(i) == 1
        expected = 1. - (1. - lmbd * 0.5 * dt) ** round(time / dt)
        self.assertAlmostEqual(stats[1][1] / stats[1][0], expected, delta = 0.08, msg = "tap root: unexpected fraction of colonized nodes")
        self.assertAlmostEqual(stats[2][1] / stats[2][0], expected, delta = 0.08, msg = "laterals: unexpected fraction of colonized nodes")


if __name__ == '__main__':
    unittest.main()
