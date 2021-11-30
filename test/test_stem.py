import unittest
import sys; sys.path.append(".."); sys.path.append("../src/python_modules")
import numpy as np
import plantbox as pb
import math


def stemAge(l, r, k, delayNG, lb):  # stem age at a certain length
    if(l == k):
        l = (1 - 1e-11) * k
    if (l < lb) or (delayNG == 0):
        return -np.log(1 - l / k) * k / r
    elif l > lb:
        return -np.log(1 - l / k) * k / r + delayNG
    else:
        raise Exception("stem length == basal length and delayNG > 0. Cannot compute exact age. ")


def stemLength(t, r, k, delayNG, lb):  # stem length at a certain age
    ageAtLb = stemAge(lb, r, k, 0, lb)  # age of stem when reach end of basal zone
    if t <= ageAtLb:
        return k * (1 - np.exp(-r * t / k))
    else:
        delay = min(delayNG, t - ageAtLb)
        return k * (1 - np.exp(-r * (t - delay) / k))


def stemLateralLength(t, et, r, k, delayNGParent, lbParent):  # length of first order laterals (without second order laterals)
    i, l = 0, 0
    if et < t:
        age = t - et
        l += stemLength(age, r, k, delayNGParent, lbParent)
    return l


class TestStem(unittest.TestCase):

    def stem_example_rtp(self, phytomereGrowth = "sequential"):
        """ an example used in the tests below, a main stem with laterals """
        self.plant = pb.Plant()  # store organism (not owned by Organ, or OrganRandomParameter)
        p0 = pb.StemRandomParameter(self.plant)
        p0.name, p0.subType, p0.la, p0.lb, p0.lmax, p0.ln, p0.r, p0.dx, p0.dxMin = "main", 1, 10., 10., 100., 1., 1.5, 1, 0.5
        p0.delayLat = 1.
        p0.delayNG = 2.
        p0.successor = [5]
        p0.successorP = [1.]

        if phytomereGrowth == "sequential":
            p0.nodalGrowth = 0
        if phytomereGrowth == "equal":
            p0.nodalGrowth = 1

        p1 = pb.StemRandomParameter(self.plant)
        p1.name, p1.subType, p1.lmax, p1.r, p1.dx, p1.dxMin = "lateral", 5, 5., 2., 1, 0.5
        self.p0, self.p1 = p0, p1  # needed at later point
        self.plant.setOrganRandomParameter(p0)  # the organism manages the type parameters and takes ownership
        self.plant.setOrganRandomParameter(p1)

        srp = pb.SeedRandomParameter(self.plant)
        self.plant.setOrganRandomParameter(srp)

        # creates seed organ (otherwise throws error in plant::simulate())
        # test == True => no need to give root parameter
        self.plant.initialize(verbose = False, test = True)

        param0 = p0.realize()  # set up stem by hand (without a stem system)
        param0.la, param0.lb = 0, 0  # its important parent has zero length, otherwise creation times are messed up
        self.ons = pb.Matrix3d(pb.Vector3d(0., 0., 1.), pb.Vector3d(0., 1., 0.), pb.Vector3d(1., 0., 0.))
        parentstem = pb.Stem(1, param0, True, True, 0., 0., self.ons, 0, False, 0)  # takes ownership of param0
        parentstem.setOrganism(self.plant)
        parentstem.addNode(pb.Vector3d(0, 0, -3), 0)  # there is no nullptr in Python
        self.parentstem = parentstem  # store parent (not owned by child Organ)
        self.stem = pb.Stem(self.plant, p0.subType, self.ons, 0, self.parentstem , 0)
        self.stem.setOrganism(self.plant)

    def stem_length_test(self, dt, l):
        ''' simulates a single stem and checks length against analytic length '''
        nl, nl2, nlth, non = [], [], [], []
        for t in dt:
            self.stem.simulate(t , False)
            nl.append(self.stem.getParameter("length"))
            nlth.append(self.stem.getParameter("lengthTh"))
            non.append(self.stem.getNumberOfNodes())

            # length from geometry
            poly = np.zeros((non[-1], 3))  #
            for i in range(0, non[-1]):
                v = self.stem.getNode(i)
                poly[i, 0] = v.x
                poly[i, 1] = v.y
                poly[i, 2] = v.z

            # relative coordinates
            sd = np.sqrt((poly ** 2).sum(axis = 1))
            nl2.append(sum(sd))
        for i in range(0, len(dt)):
            self.assertAlmostEqual(l[i], nlth[i], 10, "simulated theoretic length and analytic length do not agree in time step " + \
            str(i + 1) + " " + str(dt[i]) + " " + str(sum(dt[:(i + 1)])) + " " + str(nlth[i]) + " " + str(nl[i]) + " " + str(nl2[i]))
            self.assertAlmostEqual(nl[i], nl2[i], 10, "the two realised lengths do not agree in time step " + str(i + 1))
            self.assertLessEqual(nlth[i] - nl[i], self.p0.dxMin, "epsilonDx is too large in time step " + str(i + 1))
            self.assertLessEqual(max(sd), self.p0.dx, "axial resolution dx is too large")
            # do not take into account first node as relative coordinates == [0,0,0]
            self.assertLessEqual(self.p0.dxMin, min(sd[1:]), "axial resolution dx is unexpected small")

    def test_dynamics(self):
        ''' tests if nodes created in last time step are correct '''
        self.stem_example_rtp()
        r = self.stem
        r.simulate(2, False)
        # because of nodal growth, hasMoved() == True for leaves and stems
        self.assertEqual(r.hasMoved(), True, "dynamics: node was expected to move, but did not")
        r.simulate(1, False)
        self.assertEqual(r.hasMoved(), True, "dynamics: node was expected to move, but did not")
        non = r.getNumberOfNodes()
        r.simulate(20, False)
        self.assertEqual(r.getOldNumberOfNodes(), non, "dynamics: wrong number of old nodes")
        dx = r.getStemRandomParameter().dx
        p = r.param()
        stemLength_ = stemLength(22, p.r, p.getK(), p.delayNG, p.lb)  # r.getParameter("length")

        basalZoneLength = min(r.param().lb, stemLength_)
        basalZoneNodes = np.ceil(basalZoneLength / dx) + 1  # initial node

        branchingZoneLength = 0.
        PhytoNodes = 0.
        if(basalZoneLength < stemLength_):
            latNum = int(r.getParameter("numberOfChildren"))
            PhytoStart = basalZoneLength
            for o in range(0, latNum):
                pni = r.getChild(o).getParameter("parentNI")
                PhytoEnd = r.getLength(int(pni))
                PhytoLength = PhytoEnd - PhytoStart
                branchingZoneLength += PhytoLength
                PhytoNodes += np.ceil(PhytoLength / dx)
                PhytoStart = PhytoEnd

        apicalZoneLength = 0
        apicalZoneNodes = 0
        if((basalZoneLength + branchingZoneLength) < stemLength_):
            apicalZoneLength = stemLength_ - (basalZoneLength + branchingZoneLength)
            apicalZoneNodes = np.ceil(apicalZoneLength / dx)
        # recalculate the theoretic non
        non_th = basalZoneNodes + PhytoNodes + apicalZoneNodes
        self.assertEqual(r.getNumberOfNodes(), non_th, "dynamics: unexpected number of new nodes")  # initially, close to linear growth

    def test_leafgrow(self):
        ''' tests if the stem can create leaf '''  #
        self.stem_example_rtp()
        r = self.stem

    def test_constructors(self):
        ''' tests two kinds of constructors and copy'''
        self.stem_example_rtp()
        # 1. constructor from scratch
        param = self.p0.realize()
        stem = pb.Stem(1, param, True, True, 0., 0., self.ons , 0, False, 0)
        stem.setOrganism(self.plant)
        stem.addNode(pb.Vector3d(0, 0, -3), 0)  # parent must have at least one nodes
        # 2. used in simulation (must have parent, since there is no nullptr in Pyhton)
        stem2 = pb.Stem(self.plant, self.p1.subType, self.ons , 0, stem, 0)
        stem.addChild(stem2)
        # 3. deep copy (with a factory function)
        plant2 = pb.Organism()
        stem3 = stem.copy(plant2)
        self.assertEqual(str(stem), str(stem3), "deep copy: the organs shold be equal")
        self.assertIsNot(stem.getParam(), stem3.getParam(), "deep copy: organs have same parameter set")
        # TODO check if OTP were copied

    def test_stem_length(self):
        ''' tests if numerical stem length agrees with analytic solutions at 4 points in time'''
        self.stem_example_rtp()
        times = np.array([0., 7., 15.])  # , 30., 60.])
        dt = np.diff(times)
        k = self.stem.param().getK()  # maximal stem length
        self.assertAlmostEqual(k, 100, 12, "example stem has wrong maximal length")
        l = [stemLength(t, self.p0.r, k, self.p0.delayNG, self.p0.lb) for t in times[1:]]  # analytical stem length
        stem = self.stem.copy(self.plant)
        self.stem_length_test(dt, l)

    def test_stem_length_including_laterals(self):
        ''' tests if numerical stem length agrees with analytic solution including laterals '''
        self.stem_example_rtp()
        times = np.array([ 2., 1., 10., 15., 30., 60.])
        rp = self.stem.getStemRandomParameter()
        p = self.stem.param()  # rename
        k = p.getK()

        l = 0
        et_init = stemAge(p.lb, p.r, k, 0, p.lb)
        et = [li * p.delayLat + et_init for li in range(p.nob())]
        r2 = self.p1.r
        k2 = self.p1.lmax  # consists of lateral zone only
        t = 0.
        for i in range(0, len(times)):
            t += times[i]
            self.stem_length_test([times[i]], [stemLength(t, p.r, k, p.delayNG, p.lb)])
            organs = self.stem.getOrgans()
            nl = []
            nlth = []
            l = stemLength(t, p.r, k, p.delayNG, p.lb)
            l1 = [stemLength(max(0, t - eti), r2, k2, 0, k2) for eti in et ]
            analytic_total = l + sum(l1)
            for o in organs:
                nlth.append(o.getParameter("lengthTh"))
                nl.append(o.getParameter("length"))
            numeric_total = sum(nlth)
            self.assertAlmostEqual(numeric_total, analytic_total, 1,
            "numeric and analytic total lengths do not agree in time step " + str(i) + \
            "\nnumeric, realized: " + str(nl) + " theoretic: " + str(nlth) + " nob: " + str(len(nl) - 1) + "\nanalytice, main stem: " + str(l) + " lat: " + str(l1) + " nob: " + str(p.nob()) + \
            "\n" + str(analytic_total) + " " + str(numeric_total))

#    def test_geometry(self):
 #       """ tests if nodes can be retrieved from the organ """
        # TODO make plot for plausibility

    def test_parameter(self):
        ''' tests some parameters on sequential organ list '''
        self.stem_example_rtp()
        simtime = 30.
        self.stem.simulate(simtime, False)
        organs = self.stem.getOrgans()
        type, age, radius, order, ct = [], [], [], [], []
        for o in organs:
            type.append(o.getParameter("subType"))
            age.append(o.getParameter("age"))
            ct.append(o.getParameter("creationTime"))
            radius.append(o.getParameter("radius"))
            order.append(o.getParameter("order"))
        nol = round(self.stem.getParameter("numberOfLaterals"))
        type_ = [1.]
        type_.extend([5.] * (nol))
        self.assertEqual(type, type_, "getParameter: unexpected stem sub types")
        for i in range(0, nol):
            self.assertAlmostEqual(age[i], simtime - ct[i], 10,
            "getParameter: age and creation time does not agree for organ n°" + str(i) + \
            "\n" + str(age) + "\n" + str(ct) + "\n" + str([age[i] + ct[i] for i in range(len(ct))]))

    def test_phytomere_growth(self):
        '''check sequetial and equal growth dynamic of the phytomeres'''
        self.stem_example_rtp("equal")
        simtime = 50.
        self.stem.simulate(simtime, False)
        r = self.stem
        p = r.param()
        stemLength_th = stemLength(simtime, p.r, p.getK(), p.delayNG, p.lb)
        numPhytomeres_th = p.nob() - 1
        basalZoneLength_th = min(r.param().lb, stemLength_th)
        branchingZoneLength_th = stemLength_th - basalZoneLength_th
        PhytoLengths_th = np.full(numPhytomeres_th, branchingZoneLength_th / numPhytomeres_th)
        PhytoLengths_real = []
        stemLength_real = r.getParameter("lengthTh")
        basalZoneLength_real = r.getLength(int(r.getChild(0).getParameter("parentNI")))
        if(basalZoneLength_real < stemLength_real):
            latNum = int(r.getParameter("numberOfChildren"))
            PhytoStart = basalZoneLength_real
            for o in range(1, latNum):
                pni = r.getChild(o).getParameter("parentNI")
                PhytoEnd = r.getLength(int(pni))
                PhytoLength = PhytoEnd - PhytoStart
                PhytoLengths_real.append(PhytoLength)
                PhytoStart = PhytoEnd
        for i in range(0, len(PhytoLengths_real)):
            self.assertAlmostEqual(PhytoLengths_th[i], PhytoLengths_real[i], 10, "phytomereGrowth: unexpected phytomere length")

        self.stem_example_rtp("sequential")
        self.stem.simulate(simtime, False)
        r = self.stem
        p = r.param()
        stemLength_th = stemLength(simtime, p.r, p.getK(), p.delayNG, p.lb)
        numPhytomeres_th = p.nob() - 1
        basalZoneLength_th = min(r.param().lb, stemLength_th)
        branchingZoneLength_th = stemLength_th - basalZoneLength_th
        PhytoLengths_th = []
        PhytoLengthMin = self.p0.dxMin
        LengthToDivid = branchingZoneLength_th - PhytoLengthMin * numPhytomeres_th
        for i in range(0, numPhytomeres_th):
            PhytoLength = min(p.ln[i] - PhytoLengthMin, LengthToDivid)
            LengthToDivid -= PhytoLength
            PhytoLengths_th.append(PhytoLength + PhytoLengthMin)

        PhytoLengths_real = []
        stemLength_real = r.getParameter("lengthTh")
        basalZoneLength_real = r.getLength(int(r.getChild(0).getParameter("parentNI")))
        if(basalZoneLength_real < stemLength_real):
            latNum = int(r.getParameter("numberOfChildren"))
            PhytoStart = basalZoneLength_real
            for o in range(1, latNum):
                pni = r.getChild(o).getParameter("parentNI")
                PhytoEnd = r.getLength(int(pni))
                PhytoLength = PhytoEnd - PhytoStart
                PhytoLengths_real.append(PhytoLength)
                PhytoStart = PhytoEnd
        for i in range(0, len(PhytoLengths_real)):
            self.assertAlmostEqual(PhytoLengths_th[i], PhytoLengths_real[i], 10, "phytomereGrowth: unexpected phytomere length")


if __name__ == '__main__':
    unittest.main()

