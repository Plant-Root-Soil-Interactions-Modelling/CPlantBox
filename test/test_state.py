"""Tests for Organ::getState() / Organ::setState() on Organ and all its specialisations.

Run from this directory:
    python -m unittest test_state -v
"""

import sys

sys.path.insert(0, "..")
sys.path.insert(0, "../src/")
import json
import unittest

import plantbox as pb
from test_grassleaf import _build_plant

BASE_KEYS = {"id", "organType", "subType", "parentId", "childIds", "parentNI", "partialIHeading", "alive", "active", "age", "length", "epsilonDx", "created_linking_node", "nodes", "nodeIds", "nodeCTs", "moved", "oldNumberOfNodes"}


def _maize(days=15):
    plant = pb.Plant()
    plant.setSeed(1)
    plant.readParameters("../modelparameter/structural/plant/2020-maize.xml")
    plant.initialize(False)
    plant.simulate(days, False)
    return plant


def _grass(days=8):
    plant, _ = _build_plant()
    plant.initialize(False)
    plant.simulate(days, False)
    return plant


class TestState(unittest.TestCase):
    def organs(self, plant):
        orgs = plant.getOrgans(-1, True)
        return orgs, {o.getId(): o for o in orgs}

    def test_base_keys_and_json(self):
        plant = _maize()
        orgs, _ = self.organs(plant)
        for o in orgs:
            s = o.getState()
            self.assertTrue(BASE_KEYS <= set(s.keys()), o.getId())
            self.assertEqual(s["id"], [float(o.getId())])
            self.assertEqual(s["organType"], [float(o.organType())])
            self.assertEqual(len(s["nodes"]), 3 * o.getNumberOfNodes())
            json.loads(json.dumps(s))  # json serialisable

    def test_tree_links_are_ids(self):
        plant = _maize()
        orgs, ids = self.organs(plant)
        for o in orgs:
            s = o.getState()
            p = o.getParent()
            self.assertEqual(s["parentId"], [float(p.getId()) if p else -1.0])
            self.assertEqual(s["childIds"], [float(o.getChild(i).getId()) for i in range(o.getNumberOfChildren())])
            self.assertNotIn("plant", s)

    def test_type_specific_keys(self):
        plant = _maize()
        orgs, _ = self.organs(plant)
        orgs = orgs + [orgs[0].getParent()]  # the seed is not returned by getOrgans
        expected = {pb.seed: {"tapType", "basalType", "shootborneType", "mainStemType", "tillerType", "rootCrowns"}, pb.root: {"insertionAngle"}, pb.stem: {"stemphytomerId", "shootborneType", "localId_linking_nodes"}, pb.leaf: {"beta", "ageDependentTropism"}}
        found = set()
        for o in orgs:
            if o.organType() in expected:
                if o.organType() == pb.seed:
                    o = o.getParent() if o.getParent() else o
                found.add(o.organType())
                self.assertTrue(expected[o.organType()] <= set(o.getState().keys()))
        self.assertEqual(found, set(expected.keys()))  # seed, root, stem and leaf were all present

    def test_grassleaf_keys(self):
        plant = _grass()
        grass = [o for o in plant.getOrgans(pb.leaf, True)]
        self.assertTrue(len(grass) > 0)
        for o in grass:
            s = o.getState()
            for k in ("sheathLength", "bladeLength", "initialSheathAngle", "turtleAnchor", "turtleAnchorFrame", "turtleNodes", "turtleInitialNodeIndex"):
                self.assertIn(k, s)
            self.assertEqual(len(s["turtleNodes"]) % 4, 0)
            self.assertEqual(len(s["turtleAnchorFrame"]), 9)

    def check_restore(self, make):
        """Records the state, advances the simulation, restores the recorded state, and compares."""
        plant = make()
        orgs, ids = self.organs(plant)
        before = {o.getId(): o.getState() for o in orgs}
        plant.simulate(3, False)
        changed = sum(1 for o in orgs if o.getState() != before[o.getId()])
        self.assertGreater(changed, 0)  # the simulation changed something
        for o in orgs:
            o.setState(before[o.getId()])
        for o in orgs:
            self.assertEqual(o.getState()["age"], before[o.getId()]["age"])
            self.assertEqual(o.getState()["nodes"], before[o.getId()]["nodes"])
            self.assertEqual(o.getState()["nodeCTs"], before[o.getId()]["nodeCTs"])
            for k in before[o.getId()]:
                if k not in ("childIds",):  # laterals created later are not removed without an organ map
                    self.assertEqual(o.getState()[k], before[o.getId()][k], k)

    def test_restore_maize(self):
        self.check_restore(_maize)

    def test_restore_grassleaf(self):
        self.check_restore(_grass)

    def test_roundtrip_with_organ_map(self):
        plant = _maize()
        orgs, ids = self.organs(plant)
        for o in orgs:
            s = o.getState()
            o.setState(s, ids)
            self.assertEqual(o.getState(), s)

    def test_setstate_relinks(self):
        plant = _maize()
        orgs, ids = self.organs(plant)
        o = next(o for o in orgs if o.getNumberOfChildren() > 0)
        s = o.getState()
        s["childIds"] = s["childIds"][:1]
        o.setState(s, ids)
        self.assertEqual(o.getNumberOfChildren(), 1)
        self.assertEqual(o.getChild(0).getId(), int(s["childIds"][0]))

    def test_param_values(self):
        plant = _maize()
        orgs, _ = self.organs(plant)
        orgs = orgs + [orgs[0].getParent()]  # seed
        for o in orgs:
            s = o.getState()
            p = o.param()
            self.assertEqual(s["param.a"], [p.a])
            if o.organType() == pb.root:
                self.assertEqual(s["param.lb"], [p.lb])
                self.assertEqual(s["param.ln"], list(p.ln))
                self.assertEqual(s["param.theta"], [p.theta])
            elif o.organType() == pb.stem:
                self.assertEqual(s["param.r"], [p.r])
                self.assertEqual(s["param.nodalGrowth"], [float(p.nodalGrowth)])
            elif o.organType() == pb.leaf:
                self.assertEqual(s["param.areaMax"], [p.areaMax])
            elif o.organType() == pb.seed:
                self.assertEqual(s["param.maxB"], [float(p.maxB)])
                self.assertEqual(s["param.seedPos"], [p.seedPos.x, p.seedPos.y, p.seedPos.z])

    def test_param_grassleaf(self):
        for o in _grass().getOrgans(pb.leaf, True):
            s, p = o.getState(), o.param()
            self.assertEqual(s["param.bladeLength"], [p.bladeLength])
            self.assertEqual(s["param.sheathLength"], [p.sheathLength])
            self.assertEqual(s["param.bladeBending"], [p.bladeBending])

    def test_param_restore(self):
        plant = _maize()
        orgs, _ = self.organs(plant)
        o = next(o for o in orgs if o.organType() == pb.root)
        s = o.getState()
        s["param.lb"] = [12.5]
        s["param.ln"] = [1.0, 2.0]
        s["param.a"] = [0.7]
        o.setState(s)
        self.assertEqual(o.param().lb, 12.5)
        self.assertEqual(list(o.param().ln), [1.0, 2.0])
        self.assertEqual(o.param().a, 0.7)
        self.assertEqual(o.getState()["param.lb"], [12.5])

    def test_setstate_partial_and_id_unchanged(self):
        plant = _maize()
        o = plant.getOrgans(pb.root, True)[0]
        old = o.getState()
        o.setState({"age": [123.0], "id": [9999.0]})
        s = o.getState()
        self.assertEqual(s["age"], [123.0])
        self.assertEqual(s["id"], old["id"])
        self.assertEqual(s["nodes"], old["nodes"])


class TestOrganismState(unittest.TestCase):
    def test_keys_and_json(self):
        plant = _maize()
        s = plant.getState()
        for k in ("plantId", "simtime", "dt", "organId", "nodeId", "oldNumberOfNodes", "oldNumberOfOrgans", "rsmlSkip", "minDx", "seed_val", "stochastic", "gen", "ND", "baseOrganIds", "leafphytomerID"):
            self.assertIn(k, s)
        self.assertEqual(s["simtime"], [15.0])
        self.assertEqual(s["organId"], [float(plant.getNumberOfOrgans() - 1)])
        self.assertEqual(s["nodeId"], [float(plant.getNumberOfNodes() - 1)])
        self.assertEqual(len(s["gen"]), 625)  # mt19937: 624 words + position
        self.assertEqual(len(s["leafphytomerID"]), 10)
        self.assertEqual(s["baseOrganIds"], [float(plant.getSeed().getId())])
        json.loads(json.dumps(s))

    def test_roundtrip(self):
        plant = _maize()
        s = plant.getState()
        plant.setState(s)
        self.assertEqual(plant.getState(), s)
        orgs = plant.getOrgans(-1, True)
        plant.setState(s, {o.getId(): o for o in [plant.getSeed()] + orgs})
        self.assertEqual(plant.getState(), s)

    def test_restore_after_simulation(self):
        plant = _maize()
        s = plant.getState()
        plant.simulate(3, False)
        self.assertNotEqual(plant.getState()["simtime"], s["simtime"])
        plant.setState(s)
        self.assertEqual(plant.getState(), s)
        self.assertEqual(plant.getSimTime(), 15.0)

    def test_rng_continues_identically(self):
        a = pb.Plant()
        a.setSeed(7)
        for _ in range(5):
            a.randn()  # also fills the normal distribution cache
        s = a.getState()
        b = pb.Plant()
        b.setState(s)
        self.assertEqual([a.rand() for _ in range(5)] + [a.randn() for _ in range(5)], [b.rand() for _ in range(5)] + [b.randn() for _ in range(5)])

    def test_partial_state(self):
        plant = _maize()
        old = plant.getState()
        plant.setState({"simtime": [99.0]})
        s = plant.getState()
        self.assertEqual(s["simtime"], [99.0])
        self.assertEqual(s["gen"], old["gen"])
        self.assertEqual(s["organId"], old["organId"])


PARAM_TYPES = {pb.seed: pb.SeedSpecificParameter, pb.root: pb.RootSpecificParameter, pb.stem: pb.StemSpecificParameter, pb.leaf: pb.LeafSpecificParameter}


def all_organs(plant):
    """Seed and every organ below it, including delayed organs that have not started growing (getOrgans omits those)."""
    result, stack = [], [plant.getSeed()]
    while stack:
        o = stack.pop()
        result.append(o)
        stack.extend(o.getChild(i) for i in range(o.getNumberOfChildren()))
    return sorted(result, key=lambda o: o.getId())


def dump_plant(plant):
    """Whole plant -> JSON string: organism state plus the state of every organ (flat, linked by ids)."""
    orgs = all_organs(plant)
    return json.dumps({"organism": plant.getState(), "organs": {str(o.getId()): o.getState() for o in orgs}})


def load_plant(text, xml, cls=pb.Plant):
    """Rebuilds a plant from dump_plant(): same random parameters, organs created from their states, then linked by ids."""
    data = json.loads(text)
    plant = cls()
    plant.readParameters(xml)
    plant.initialize(False)  # callbacks and seed-derived parameters; the created organs are replaced below
    organs = {}
    for key, s in sorted(data["organs"].items(), key=lambda kv: int(kv[0])):
        ot, oid = int(s["organType"][0]), int(s["id"][0])
        param = PARAM_TYPES[ot]()
        if ot == pb.seed:
            o = pb.Seed(oid, param, True, True, 0.0, 0.0, False, 0)
        else:
            ctor = {pb.root: pb.Root, pb.stem: pb.Stem, pb.leaf: pb.Leaf}[ot]
            o = ctor(oid, param, True, True, 0.0, 0.0, pb.Vector3d(0, 0, 1), 0, False, 0)
        o.setOrganism(plant)
        organs[oid] = o
    for oid, o in organs.items():
        o.setState(data["organs"][str(oid)], organs)
        p = o.getState()["parentId"][0]
        if p >= 0:
            o.setParent(organs[int(p)])
    plant.setState(data["organism"], organs)
    return plant


class TestWholePlant(unittest.TestCase):
    XML = "../modelparameter/structural/plant/2020-maize.xml"

    def compare(self, a, b):
        oa, ob = all_organs(a), all_organs(b)
        self.assertEqual(len(oa), len(ob))
        self.assertEqual(a.getState(), b.getState())
        for x, y in zip(oa, ob):
            sx, sy = x.getState(), y.getState()
            for k in sx:
                self.assertEqual(sx[k], sy[k], (x.getId(), k))

    def test_dump_load_identical(self):
        plant = _maize()
        copy = load_plant(dump_plant(plant), self.XML)
        self.compare(plant, copy)
        self.assertEqual(copy.getNumberOfNodes(), plant.getNumberOfNodes())
        self.assertEqual(len(copy.getSegments()), len(plant.getSegments()))

    def test_continue_simulation_identically(self):
        plant = _maize()
        copy = load_plant(dump_plant(plant), self.XML)
        plant.simulate(5, False)
        copy.simulate(5, False)
        self.assertGreater(plant.getNumberOfOrgans(), 0)
        self.compare(plant, copy)


if __name__ == "__main__":
    unittest.main()
