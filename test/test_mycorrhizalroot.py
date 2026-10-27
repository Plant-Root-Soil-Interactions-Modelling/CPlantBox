import sys; sys.path.append(".."); sys.path.append("../src/")
import unittest

import plantbox as pb

import numpy as np
from scipy.linalg import norm


def rootAge(l, r, k):  # root age at a certain length
    return -np.log(1 - l / k) * k / r


def rootLength(t, r, k):  # root length at a certain age
    return k * (1 - np.exp(-r * t / k))


def rootLateralLength(t, et, r, k):  # length of first order laterals (without second order laterals)
    i, l = 0, 0
    while et[i] < t:
        age = t - et[i]
        l += rootLength(age, r, k)
        i += 1
    return l


def segment_lengths(organ):
    """ lengths of all segments of the organ """
    nodes = np.array([[n.x, n.y, n.z] for n in organ.getNodes()])
    return norm(np.diff(nodes, axis = 0), axis = 1)


def colonization(root):
    """ colonization state and time of all nodes of the root """
    n = root.getNumberOfNodes()
    return np.array([root.getNodeColonization(i) for i in range(n)]), np.array([root.getNodeColonizationTime(i) for i in range(n)])


class TestMycorrhizalRoot(unittest.TestCase):
    """ All tests simulate a MycorrhizalPlant with a seed: MycorrhizalRoot limits node creation times to the end of the
        current time step of the plant, and hyphae need a MycorrhizalPlant, so a root cannot be simulated on its own """

    def mycroot_plant(self, lmbd = 0., vi = 0., minAge = 0., maxAge = 1.e6, highresolution = 0., dx_inf = 0.1, hed = 0.,
                      hyphalDelay = 0., laterals = True, delay_definition = None, seed = 1):
        """ a mycorrhizal tap root with laterals, without colonization (lmbd = 0) by default """
        self.plant = pb.MycorrhizalPlant(seed)
        srp = pb.SeedRandomParameter(self.plant)
        srp.subType = 0
        if delay_definition is not None:
            srp.delayDefinitionRoot = delay_definition
        self.plant.setOrganRandomParameter(srp)
        p0 = pb.MycorrhizalRootRandomParameter(self.plant)
        p0.name, p0.subType, p0.la, p0.lb, p0.lmax, p0.ln, p0.r, p0.dx = "taproot", 1, 10., 1., 100., 1., 1.5, 0.5
        if laterals:
            p0.successor = [[2]]
            p0.successorP = [[1.]]
        p1 = pb.MycorrhizalRootRandomParameter(self.plant)
        p1.name, p1.subType, p1.lmax, p1.r, p1.dx = "lateral", 2, 25., 2., 0.1
        for p in [p0, p1]:
            p.lmbd, p.vi, p.minAge, p.maxAge, p.highresolution, p.dx_inf = lmbd, vi, minAge, maxAge, highresolution, dx_inf
            p.hyphalEmergenceDensity, p.hyphalDelay = hed, hyphalDelay
            self.plant.setOrganRandomParameter(p)
        self.p0, self.p1 = p0, p1
        self.add_hyphae_parameter()
        self.plant.initialize(False)
        self.root = self.plant.getSeed().getChild(0)  # the tap root (getOrgans() skips it, as long as it has only one node)

    def fast_root_plant(self, lmbd = 0., vi = 0., minAge = 0., maxAge = 1.e6, highresolution = 0., dx_inf = 0.1, hed = 0.,
                        hyphalDelay = 0., seed = 1):
        """ a mycorrhizal tap root without laterals, growing linearly to its final length of 10 cm within the first time step of 0.1 day,
            so that colonization takes place on a root of constant geometry """
        self.plant = pb.MycorrhizalPlant(seed)
        srp = pb.SeedRandomParameter(self.plant)
        srp.subType = 0
        self.plant.setOrganRandomParameter(srp)
        p0 = pb.MycorrhizalRootRandomParameter(self.plant)
        p0.name, p0.subType, p0.lmax, p0.r, p0.dx, p0.gf = "taproot", 1, 10., 100., 0.5, 2  # linear growth
        p0.lmbd, p0.vi, p0.minAge, p0.maxAge, p0.highresolution, p0.dx_inf = lmbd, vi, minAge, maxAge, highresolution, dx_inf
        p0.hyphalEmergenceDensity, p0.hyphalDelay = hed, hyphalDelay
        self.plant.setOrganRandomParameter(p0)
        self.p0 = p0
        self.add_hyphae_parameter()
        self.plant.initialize(False)
        self.root = self.plant.getSeed().getChild(0)

    def add_hyphae_parameter(self):
        """ root hyphae are of sub type 1, short and unbranched to keep the tests fast """
        self.hrp = pb.HyphaeRandomParameter(self.plant)
        self.hrp.name, self.hrp.subType, self.hrp.dx, self.hrp.lmax, self.hrp.v = "hyphae", 1, 0.05, 0.2, 0.1
        self.plant.setOrganRandomParameter(self.hrp)

    def simulate(self, time, dt = 0.1):
        for i in range(0, int(round(time / dt))):
            self.plant.simulate(dt, False)

    def root_hyphae(self, root):
        return [root.getChild(i) for i in range(root.getNumberOfChildren()) if root.getChild(i).organType() == pb.hyphae]

    def root_laterals(self, root):
        return [root.getChild(i) for i in range(root.getNumberOfChildren()) if root.getChild(i).organType() == pb.root]

    """ root geometry (as in test_root.py) """

    def test_constructors(self):
        """ tests two kinds of constructors and copy """
        self.mycroot_plant()
        param = self.p0.realize()
        root = pb.MycorrhizalRoot(1, param, True, True, 0., 0., pb.Vector3d(0, 0, -1), 0, False, 0)  # 1. constructor from scratch
        root.setOrganism(self.plant)
        root.addNode(pb.Vector3d(0, 0, -3), 0)  # parent must have at least one nodes
        root2 = pb.MycorrhizalRoot(self.plant, self.p1.subType, 0, root, 0)  # 2. used in simulation (must have parent)
        root.addChild(root2)
        plant2 = pb.MycorrhizalPlant()
        root3 = root.copy(plant2)  # 3. deep copy
        self.assertEqual(root3.__class__.__name__, "MycorrhizalRoot", "deep copy: unexpected class type")
        self.assertEqual(str(root), str(root3), "deep copy: the root string representations shold be equal")
        self.assertIsNot(root.getParam(), root3.getParam(), "deep copy: roots have same specific parameter set")
        self.assertEqual(str(root.param()), str(root3.param()), "deep copy: roots have different parameter values")

    def test_root_length(self):
        """ tests if numerical root length agrees with analytic solutions at 4 points in time with two scales of dt """
        times = np.array([0., 7., 15., 30., 60.])
        for subDt in [1, 100]:
            self.mycroot_plant(laterals = False)
            k = self.root.param().getK()
            self.assertAlmostEqual(k, 100, 12, "example root has wrong maximal length")
            for i, dt in enumerate(np.diff(times)):
                self.simulate(dt, dt / subDt)
                l = rootLength(times[i + 1], self.p0.r, k)
                self.assertAlmostEqual(self.root.getParameter("length"), l, 10, "numeric and analytic lengths do not agree in time step " + str(i + 1))
                self.assertAlmostEqual(sum(segment_lengths(self.root)), l, 10, "numeric and geometric lengths do not agree in time step " + str(i + 1))

    def test_root_length_including_laterals(self):
        """ tests if numerical root length agrees with analytic solution including laterals """
        times = np.array([0., 7., 15., 30., 60.])
        for subDt in [1, 100]:
            self.mycroot_plant()
            rp = self.root.getRootRandomParameter()
            p = self.root.param()
            k = p.getK()
            et = np.zeros((p.nob()))  # lateral emergence times
            l = 0
            et[0] = rootAge(p.la - rp.ln / 2. + p.lb + l, p.r, k)
            for i in range(0, p.nob() - 1):
                l += p.ln[i]
                et[i + 1] = rootAge(p.la - rp.ln / 2. + p.lb + l, p.r, k)
            for i, dt in enumerate(np.diff(times)):
                self.simulate(dt, dt / subDt)
                t = times[i + 1]
                analytic = rootLength(t, p.r, k) + rootLateralLength(t, et, self.p1.r, self.p1.lmax)
                numeric = sum([o.getParameter("length") for o in self.plant.getOrgans(pb.root)])
                self.assertAlmostEqual(numeric, analytic, 10, "numeric and analytic total lengths do not agree in time step " + str(i + 1))

    def test_parameter(self):
        """ tests some parameters on sequential organ list """
        self.mycroot_plant()
        simtime = 30.
        self.simulate(simtime, 1.)
        organs = self.plant.getOrgans(pb.root)
        nol = round(self.root.getParameter("numberOfLaterals"))
        self.assertEqual([o.getParameter("subType") for o in organs], [1.] + [2.] * nol, "getParameter: unexpected root sub types")
        order0 = self.root.getParameter("order")
        for o in organs[1:]:
            self.assertEqual(o.getParameter("order"), order0 + 1, "getParameter: unexpected root order")
        for o in organs:
            self.assertAlmostEqual(o.getParameter("age"), simtime - o.getParameter("creationTime"), 10, "getParameter: age and creation time does not agree")

    def test_dynamics(self):
        """ tests if nodes created in last time step are correct """
        self.mycroot_plant()
        r = self.root
        self.plant.simulate(.5, False)
        self.assertEqual(r.hasMoved(), False, "dynamics: node is created during first step")
        self.plant.simulate(1e-1, False)
        self.assertEqual(r.hasMoved(), True, "dynamics: node was expected to move, but did not")
        non = r.getNumberOfNodes()
        self.plant.simulate(2.4, False)
        self.assertEqual(r.getOldNumberOfNodes(), non, "dynamics: wrong number of old nodes")
        self.assertEqual(r.getNumberOfNodes() - non, round(2.4 * r.param().r / self.p0.dx), "dynamics: unexpected number of new nodes")  # initially, close to linear growth

    def test_new_delay_types(self):
        """ tests the emergence time of laterals for the three delay definitions """
        self.mycroot_plant(delay_definition = 0)  # depends on la
        time = 100
        self.plant.simulate(time, False)
        effectiveLa = max(self.root.getParameter("la") - self.root.getParameter("lnMean") / 2, 0.)
        for i, r in enumerate(self.root_laterals(self.root)):
            lengthAtCreation = self.root.getLength(r.parentNI)
            et = self.root.calcAge(lengthAtCreation + effectiveLa)  # age of the root, when the apical zone is developed
            self.assertAlmostEqual(r.getAge(), time - et, 10, "numeric and analytic age of root n#" + str(i + 1) + " do not agree")
        for delay_definition in [1, 2]:  # depends on ldelay of the parent organ (1) or of the lateral (2)
            self.mycroot_plant(delay_definition = delay_definition)
            self.p0.ldelay, self.p0.ldelays = 1., 0.
            self.p1.ldelay, self.p1.ldelays = 2., 0.
            time = 10
            self.plant.simulate(time, False)
            for i, r in enumerate(self.root_laterals(self.root)):
                dl = self.p0.ldelay if delay_definition == 1 else self.p1.ldelay
                et = self.root.getNodeCT(r.parentNI) + dl
                self.assertAlmostEqual(r.getAge(), time - et, 10, "numeric and analytic age of root n#" + str(i + 1) + " do not agree")

    """ colonization """

    def test_no_colonization(self):
        """ tests that nothing is colonized for lmbd = 0 """
        self.mycroot_plant(vi = 0.5, hed = 1.)
        self.simulate(10., 1.)
        for r in self.plant.getOrgans(pb.root):
            col, _ = colonization(r)
            self.assertTrue(np.all(col == 0), "no colonization: root is colonized although lmbd = 0")
        self.assertEqual(len(self.plant.getOrgans(pb.hyphae)), 0, "no colonization: hyphae emerged although nothing is colonized")

    def test_primary_colonization(self):
        """ tests the mean fraction of primarily colonized nodes against the colonization probability lmbd*segment length*dt per time step """
        lmbd, dt, time = 0.2, 0.1, 5.
        n, n_col = 0, 0
        for seed in range(1, 11):
            self.fast_root_plant(lmbd = lmbd, seed = seed)
            self.simulate(time, dt)
            col, _ = colonization(self.root)
            self.assertFalse(np.any(col > 1), "primary colonization: secondary colonization although vi = 0")
            n += len(col) - 1  # node 0 cannot be colonized primarily
            n_col += np.sum(col[1:] == 1)
            seg = segment_lengths(self.root)
            self.assertAlmostEqual(self.root.getParameter("primaryColonization"), np.sum(seg[col[1:] == 1]), 10,
                                   "primary colonization: getParameter('primaryColonization') unexpected")
        p_expected = 1. - (1. - lmbd * self.p0.dx * dt) ** round(time / dt)
        self.assertAlmostEqual(n_col / n, p_expected, delta = 0.1, msg = "primary colonization: unexpected fraction of colonized nodes")

    def test_colonization_age(self):
        """ tests that nodes are only colonized primarily at an age between minAge and maxAge """
        minAge, maxAge = 1., 2.
        self.fast_root_plant(lmbd = 5., minAge = minAge, maxAge = maxAge)
        self.simulate(4.)
        col, ct = colonization(self.root)
        self.assertGreater(np.sum(col == 1), 0, "colonization age: no colonization")
        for i in np.where(col == 1)[0]:
            age = ct[i] - self.root.getNodeCT(i)
            self.assertGreaterEqual(age, minAge - 1.e-9, "colonization age: node colonized before minAge")
            self.assertLessEqual(age, maxAge + 1.e-9, "colonization age: node colonized after maxAge")

    def test_secondary_colonization(self):
        """ tests that secondary colonization spreads from colonized nodes with speed vi at most """
        vi = 0.5
        self.fast_root_plant(lmbd = 0.05, vi = vi, seed = 3)
        self.simulate(5.)
        col, ct = colonization(self.root)
        seg = segment_lengths(self.root)
        self.assertGreater(np.sum(col == 2), 0, "secondary colonization: no secondary colonization")
        for j in np.where(col == 2)[0]:
            neighbours = [k for k in [j - 1, j + 1] if 0 <= k < len(col) and col[k] > 0]
            self.assertGreater(len(neighbours), 0, "secondary colonization: isolated secondary colonization at node " + str(j))
            earliest = min([ct[k] + seg[min(j, k)] / vi for k in neighbours])  # the front needs at least seg/vi from a neighbour
            self.assertGreaterEqual(ct[j], earliest - 1.e-9, "secondary colonization: colonization spread faster than vi at node " + str(j))
            self.assertLessEqual(ct[j], self.plant.getSimTime() + 1.e-9, "secondary colonization: colonization time in the future")
        self.assertAlmostEqual(self.root.getParameter("secondaryColonization"), np.sum(seg[col[1:] > 1]), 10,
                               "secondary colonization: getParameter('secondaryColonization') unexpected")
        self.assertAlmostEqual(self.root.getParameter("colonizationLength"), np.sum(seg[col[1:] > 0]), 10,
                               "secondary colonization: getParameter('colonizationLength') unexpected")

    def test_colonization_transfer(self):
        """ tests that colonization is transferred from the parent root to the base of its laterals (state 3) """
        self.mycroot_plant(lmbd = 0.1, vi = 0.5)
        self.simulate(20., 0.5)
        col, ct = colonization(self.root)
        transfers = 0
        for l in self.root_laterals(self.root):
            if l.getNumberOfNodes() < 2:
                continue
            c0 = l.getNodeColonization(0)
            if c0 == 3:
                transfers += 1
                self.assertGreater(col[l.parentNI], 0, "colonization transfer: parent node of a transferred lateral is not colonized")
                self.assertGreaterEqual(l.getNodeColonizationTime(0), max(l.getNodeCT(0), ct[l.parentNI]) - 1.e-9,
                                        "colonization transfer: lateral colonized before its parent node or its creation")
            elif col[l.parentNI] > 0:  # transfer is skipped only if the first segment of the lateral is already colonized
                self.assertGreater(l.getNodeColonization(1), 0, "colonization transfer: colonization was not transferred to the lateral")
        self.assertGreater(transfers, 0, "colonization transfer: no colonization transfer")

    def test_colonization_persistent(self):
        """ tests that colonization state and time of a node do not change once it is colonized (also with node refinement);
            nodes are identified by their position, inserted nodes get new global node indices, as in Stem::addNode """
        self.fast_root_plant(lmbd = 0.2, vi = 0.2, highresolution = 1.)
        self.simulate(1.)

        def by_position():
            col, ct = colonization(self.root)
            return {tuple(np.round([n.x, n.y, n.z], 9)): (col[i], ct[i]) for i, n in enumerate(self.root.getNodes())}

        before = {pos: v for pos, v in by_position().items() if v[0] > 0}
        self.assertGreater(len(before), 0, "colonization persistent: no colonization")
        self.simulate(2.)
        after = by_position()
        for nid, (c, t) in before.items():
            self.assertIn(nid, after, "colonization persistent: colonized node disappeared")
            self.assertEqual(after[nid][0], c, "colonization persistent: colonization state changed")
            self.assertEqual(after[nid][1], t, "colonization persistent: colonization time changed")

    @unittest.expectedFailure
    def test_segment_refinement(self):
        """ tests that colonized segments are refined to dx_inf (highresolution = 1) without changing the root length;
            fails: MycorrhizalRoot::addNode inserts the new node at index-1, i.e. before node i-1 instead of between i-1 and i """
        dx_inf = 0.1
        self.fast_root_plant(lmbd = 0.2, vi = 0.2, highresolution = 1., dx_inf = dx_inf)
        self.simulate(3.)
        col, ct = colonization(self.root)
        seg = segment_lengths(self.root)
        self.assertEqual(len(col), self.root.getNumberOfNodes(), "segment refinement: colonization vector has wrong size")
        self.assertGreater(self.root.getNumberOfNodes(), self.p0.lmax / self.p0.dx + 1, "segment refinement: no nodes were inserted")
        self.assertAlmostEqual(self.root.getLength(), self.p0.lmax, 10, "segment refinement: root length changed")
        self.assertAlmostEqual(np.sum(seg), self.p0.lmax, 10, "segment refinement: geometric root length changed")
        for i in range(1, len(col)):
            if col[i] > 0 and col[i - 1] > 0:
                self.assertLessEqual(seg[i - 1], dx_inf + 1.e-9, "segment refinement: colonized segment longer than dx_inf")
        for h in self.root_hyphae(self.root):  # children have to be moved with the inserted nodes
            self.assertEqual(h.getNodeId(0), self.root.getNodeId(h.parentNI), "segment refinement: hypha is not attached to its root node")

    """ hyphal emergence """

    def test_hyphal_emergence_density(self):
        """ tests that hyphae emerge every 1/hyphalEmergenceDensity cm along the colonized root (highresolution = 0) """
        hed = 2.
        self.fast_root_plant(lmbd = 1000., hed = hed)  # root is fully colonized in the first time step
        self.simulate(1.)
        col, ct = colonization(self.root)
        self.assertTrue(np.all(col[1:] > 0), "hyphal emergence: root should be fully colonized")
        hyphae = self.root_hyphae(self.root)
        self.assertEqual(len(hyphae), int(np.floor(hed * self.root.getLength() + 1.e-9)), "hyphal emergence: unexpected number of hyphae")
        for h in hyphae:
            self.assertGreater(col[h.parentNI], 0, "hyphal emergence: hypha emerged from a node that is not colonized")
            self.assertEqual(h.getNodeId(0), self.root.getNodeId(h.parentNI), "hyphal emergence: hypha is not attached to its root node")

    def test_hyphal_emergence_highresolution(self):
        """ tests that exactly one hypha emerges at every colonized node (highresolution = 1) """
        self.fast_root_plant(lmbd = 0.2, vi = 0.2, highresolution = 1.)
        self.simulate(2.)
        col, ct = colonization(self.root)
        pnis = [h.parentNI for h in self.root_hyphae(self.root)]
        self.assertGreater(len(pnis), 0, "hyphal emergence: no hyphae")
        for i in range(1, len(col)):
            self.assertEqual(pnis.count(i), 1 if col[i] > 0 else 0, "hyphal emergence: unexpected number of hyphae at node " + str(i))

    def test_hyphal_delay(self):
        """ tests that hyphae emerge hyphalDelay after the colonization of their root node """
        hyphalDelay = 0.5
        self.fast_root_plant(lmbd = 0.5, hed = 2., hyphalDelay = hyphalDelay)
        self.simulate(3.)
        col, ct = colonization(self.root)
        hyphae = self.root_hyphae(self.root)
        self.assertGreater(len(hyphae), 0, "hyphal delay: no hyphae")
        for h in hyphae:
            self.assertAlmostEqual(h.getNodeCT(0), ct[h.parentNI] + hyphalDelay, 9, "hyphal delay: unexpected emergence time")
            self.assertLessEqual(h.getAge(), self.plant.getSimTime() - (ct[h.parentNI] + hyphalDelay) + 1.e-9, "hyphal delay: hypha too old")

    """ other """

    def test_set_colonization(self):
        """ tests setting the colonization of a node manually """
        self.fast_root_plant()
        self.simulate(0.1)
        self.root.setColonization(3, 1, 0.05)
        self.assertEqual(self.root.getNodeColonization(3), 1, "setColonization: unexpected state")
        self.assertEqual(self.root.getNodeColonizationTime(3), 0.05, "setColonization: unexpected time")
        with self.assertRaises(IndexError):
            self.root.getNodeColonization(self.root.getNumberOfNodes())

    def test_copy_colonization(self):
        """ tests that a deep copy keeps the colonization of the nodes """
        self.fast_root_plant(lmbd = 0.2, vi = 0.2, hed = 1.)
        self.simulate(2.)
        root2 = self.root.copy(pb.MycorrhizalPlant())
        col, ct = colonization(self.root)
        col2, ct2 = colonization(root2)
        self.assertTrue(np.array_equal(col, col2), "copy: colonization states differ")
        self.assertTrue(np.array_equal(ct, ct2), "copy: colonization times differ")
        self.assertEqual(len(self.root_hyphae(root2)), len(self.root_hyphae(self.root)), "copy: hyphae are not copied")


if __name__ == '__main__':
    unittest.main()
