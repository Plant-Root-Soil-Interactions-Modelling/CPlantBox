import unittest
import sys
import numpy as np
sys.path.append("..")
import plantbox as pb


def stemAge(l, r, k):  # stem age at a certain length
    return -np.log(1 - l / k) * k / r


def stemLength(t, r, k):  # stem length at a certain age
    return k * (1 - np.exp(-r * t / k))


def stemLateralLength(t, et, r, k):  # length of first order laterals (without second order laterals)
    i, l = 0, 0
    while et[i] < t:
        age = t - et[i]
        l += stemLength(age, r, k)
        i += 1
    return l


class TestStem(unittest.TestCase):

    def stem_example_rtp(self):
        """ an example used in the tests below, a main stem with laterals """
        self.plant = pb.Organism()  # store organism (not owned by Organ, or OrganRandomParameter)
        p0 = pb.StemRandomParameter(self.plant)
        p0.name, p0.subType, p0.la, p0.lb, p0.lmax, p0.ln, p0.r, p0.dx = "main", 1,  10., 1., 100., 1., 1.5, 0.5
        p0.successor = [5]
        p0.successorP = [1.]
        p1 = pb.StemRandomParameter(self.plant)
        p1.name, p1.subType, p1.lmax, p1.r, p1.dx = "lateral", 5, 25., 2., 0.1
        self.p0, self.p1 = p0, p1  # needed at later point
        self.plant.setOrganRandomParameter(p0)  # the organism manages the type parameters and takes ownership
        self.plant.setOrganRandomParameter(p1)
        # TODO (first node is not set, if seed is used)
        srp = pb.SeedRandomParameter(self.plant)
        self.plant.setOrganRandomParameter(srp)
        #
        param0 = p0.realize()  # set up stem by hand (without a stem system)
        param0.la, param0.lb = 0, 0  # its important parent has zero length, otherwise creation times are messed up
        parentstem = pb.Stem(1, param0, True, True, 0., 0., pb.Vector3d(0, 0, -1), 0, 0, False, 0)  # takes ownership of param0
        parentstem.setOrganism(self.plant)
        parentstem.addNode(pb.Vector3d(0, 0, -3), 0)  # there is no nullptr in Python
        self.parentstem = parentstem  # store parent (not owned by child Organ)
        #
        self.stem = pb.Stem(self.plant, p0.subType, pb.Vector3d(0, 0, -1), 0, self.parentstem , 0, 0)
        self.stem.setOrganism(self.plant)

    def stem_length_test(self, dt, l, subDt):
        """ simulates a single stem and checks length against analytic length """
        nl, nl2, non, meanDX = [], [], [], []
        for t in dt:
            for i in range(0, subDt):

                self.stem.simulate(t / subDt, False)
            nl.append(self.stem.getParameter("length"))
            non.append(self.stem.getNumberOfNodes())
            meanDX.append(nl[-1] / non[-1])

            # length from geometry
            poly = np.zeros((non[-1], 3))  #
            for i in range(0, non[-1]):
                v = self.stem.getNode(i)
                poly[i, 0] = v.x
                poly[i, 1] = v.y
                poly[i, 2] = v.z
            d = np.diff(poly, axis = 0)
            sd = np.sqrt((d ** 2).sum(axis = 1))
            nl2.append(sum(sd))
        for i in range(0, len(dt)):
            self.assertAlmostEqual(l[i], nl[i], 10, "numeric and analytic lengths do not agree in time step " + str(i + 1))
            self.assertAlmostEqual(l[i], nl2[i], 10, "numeric and analytic lengths do not agree in time step " + str(i + 1))
            self.assertLessEqual(meanDX[i], 0.5, "axial resolution dx is too large")
            self.assertLessEqual(0.25, meanDX[i], "axial resolution dx is unexpected small")

    def test_constructors(self):
        """ tests two kinds of constructors and copy"""
        self.stem_example_rtp()
        # 1. constructor from scratch
        param = self.p0.realize()
        stem = pb.Stem(1, param, True, True, 0., 0., pb.Vector3d(0, 0, -1), 0, 0, False, 0)
        stem.setOrganism(self.plant)
        stem.addNode(pb.Vector3d(0, 0, -3), 0)  # parent must have at least one nodes
        # 2. used in simulation (must have parent, since there is no nullptr in Pyhton)
        stem2 = pb.Stem(self.plant, self.p1.subType, pb.Vector3d(0, 0, -1), 0, stem, 0, 0)
        stem.addChild(stem2)
        # 3. deep copy (with a factory function)
        plant2 = pb.Organism()
        stem3 = stem.copy(plant2)
        self.assertEqual(str(stem), str(stem3), "deep copy: the organs shold be equal")
        self.assertIsNot(stem.getParam(), stem3.getParam(), "deep copy: organs have same parameter set")
        # TODO check if OTP were copied

    def test_stem_length(self):
        """ tests if numerical stem length agrees with analytic solutions at 4 points in time with two scales of dt"""
        self.stem_example_rtp()
        times = np.array([0., 7., 15., 30., 60.])
        dt = np.diff(times)
        k = self.stem.param().getK()  # maximal stem length
        self.assertAlmostEqual(k, 100, 12, "example stem has wrong maximal length")
        l = stemLength(times[1:], self.p0.r, k)  # analytical stem length
        stem = self.stem.copy(self.plant)
        self.stem_length_test(dt, l, 1)  # large dt
        self.stem = stem
        self.stem_length_test(dt, l, 1000)  # very fine dt

    def test_stem_length_including_laterals(self):
        """ tests if numerical stem length agrees with analytic solution including laterals """
        self.stem_example_rtp()
        times = np.array([0., 7., 15., 30., 60.])
        dt = np.diff(times)
        rp = self.stem.getStemRandomParameter()                                       
        p = self.stem.param()  # rename
        k = p.getK()
        print(k)
        et = np.zeros((p.nob()))
        l = 0
        et[0] = stemAge(p.la - rp.ln / 2. + p.lb + l, p.r, k)
        for i in range(0, p.nob() - 1):  # calculate lateral emergence times
            l += p.ln[i]
            et[i + 1] = stemAge(p.la - rp.ln / 2. + p.lb + l, p.r, k)
        l = stemLength(times[1:], p.r, k)  # zero order lengths
        l1 = []
        r2 = self.p1.r
        k2 = self.p1.lmax  # consists of lateral zone only
        for t in times[1:]:
            l1.append(stemLateralLength(t, et, r2, k2))
        analytic_total = l + l1

        for subDX in [1, 1000]:
            numeric_total = []
            for t in times[1:]:
                stem = self.stem.copy(self.plant)
                self.stem_length_test([t], [stemLength(t, p.r, k)], subDX)
                organs = self.stem.getOrgans()
                nl = 0
                for o in organs:
                    nl += o.getParameter("length")
                numeric_total.append(nl);
                self.stem = stem
            for i in range(0, len(times[1:])):
                self.assertAlmostEqual(numeric_total[i], analytic_total[i], 10, "numeric and analytic total lengths do not agree in time step " + str(i + 1))

    def test_geometry(self):
        """ tests if nodes can be retrieved from the organ """
        # TODO make plot for plausibility

    def test_parameter(self):
        """ tests some parameters on sequential organ list """
        self.stem_example_rtp()
        simtime = 30.             
        self.stem.simulate(30,False)
        organs = self.stem.getOrgans()
        type, age, radius, order, ct = [], [], [], [], []
        for o in organs:
            type.append(o.getParameter("subType"))
            age.append(o.getParameter("age"))
            ct.append(o.getParameter("creationTime"))
            radius.append(o.getParameter("radius"))
            order.append(o.getParameter("order"))
#        nol = round(self.stem.getParameter("numberOfLaterals"))
        type_ = [1.]
        type_.extend([2.] * nol)                                                       
        self.assertEqual(type, type_, "getParameter: unexpected stem sub types")
        self.assertEqual(order, type_, "getParameter: unexpected stem order")  # +1, because of artificial parent root
        for i in range(0, nol):
            self.assertAlmostEqual(age[i], simtime - ct[i], 10, "getParameter: age and creation time does not agree") 
            
                        
       
    def test_dynamics(self):
        """ tests if nodes created in last time step are correct """  #
        self.stem_example_rtp()
        r = self.stem
        r.simulate(.5, True)
        self.assertEqual(r.hasMoved(), False, "dynamics: node movement during first step")
        r.simulate(.1, True)
                          
        self.assertEqual(r.hasMoved(), True, "dynamics: node was expected to move, but did not")
        non = r.getNumberOfNodes() 
        r.simulate(2.4, False)
        self.assertEqual(r.getOldNumberOfNodes(), non, "dynamics: wrong number of old nodes")
        dx = r.getStemRandomParameter().dx
        self.assertEqual(r.getNumberOfNodes() - non, round(2.4 * r.param().r / dx), "dynamics: unexpected number of new nodes")  # initially, close to linear growth

                                 
        
    def test_leafgrow(self):
        """ tests if the stem can create leaf """  #
        self.stem_example_rtp()
        r = self.stem

if __name__ == '__main__':
    unittest.main()
