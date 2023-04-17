import sys; sys.path.append(".."); sys.path.append("../src/")
import unittest

import plantbox as pb

import numpy as npfrom scipy.linalg import norm


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


class TestRoot(unittest.TestCase):

    def root_example_rrp(self):
        """ an example used in the tests below, a main root with laterals """
        self.plant = pb.Organism()  # store organism (not owned by Organ, or OrganRandomParameter)
        p0 = pb.RootRandomParameter(self.plant)
        p0.name, p0.subType, p0.la, p0.lb, p0.lmax, p0.ln, p0.r, p0.dx = "taproot", 1, 10., 1., 100., 1., 1.5, 0.5
        p0.successor = [[2]]
        p0.successorP = [[1.]]
        p1 = pb.RootRandomParameter(self.plant)
        p1.name, p1.subType, p1.lmax, p1.r, p1.dx = "lateral", 2, 25., 2., 0.1
        self.p0, self.p1 = p0, p1  # needed at later point
        self.plant.setOrganRandomParameter(p0)  # the organism manages the type parameters and takes ownership
        self.plant.setOrganRandomParameter(p1)
        srp = pb.SeedRandomParameter(self.plant)
        self.plant.setOrganRandomParameter(srp)

        param0 = p0.realize()  # set up root by hand (without a root system)
        param0.la, param0.lb = 0, 0  # its important parent has zero length, otherwise creation times are messed up
        parentroot = pb.Root(1, param0, True, True, 0., 0., pb.Vector3d(0, 0, -1), 0, False, 0)  # takes ownership of param0
        parentroot.setOrganism(self.plant)
        parentroot.addNode(pb.Vector3d(0, 0, -3), 0)  # there is no nullptr in Python

        self.parentroot = parentroot  # store parent (not owned by child Organ)
        self.root = pb.Root(self.plant, p0.subType,  0, self.parentroot , 0)
        self.root.setOrganism(self.plant)

    def root_length_test(self, dt, l, subDt):
        """ simulates a single root and checks length against analytic length """
        nl, nl2, non, meanDX = [], [], [], []
        for t in dt:
            for i in range(0, subDt):
                self.root.simulate(t / subDt, False)
            nl.append(self.root.getParameter("length"))
            non.append(self.root.getNumberOfNodes())
            meanDX.append(nl[-1] / non[-1])
            # length from geometry
            poly = np.zeros((non[-1], 3))  #
            for i in range(0, non[-1]):
                v = self.root.getNode(i)
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
        self.root_example_rrp()
        # 1. constructor from scratch
        param = self.p0.realize()
        root = pb.Root(1, param, True, True, 0., 0., pb.Vector3d(0, 0, -1), 0, False, 0)
        root.setOrganism(self.plant)
        root.addNode(pb.Vector3d(0, 0, -3), 0)  # parent must have at least one nodes
        # 2. used in simulation (must have parent, since there is no nullptr in Pyhton)
        root2 = pb.Root(self.plant, self.p1.subType, 0, root, 0)
        root.addChild(root2)
        # 3. deep copy (with a factory function)
        plant2 = pb.Organism()
        root3 = root.copy(plant2)
        self.assertEqual(str(root), str(root3), "deep copy: the root string representations shold be equal")
        self.assertIsNot(root.getParam(), root3.getParam(), "deep copy: roots have same specific parameter set")  # type OrganSpecificParameter
        self.assertEqual(str(root.param()), str(root3.param()), "deep copy: roots have different parameter values")  # type RootSpecificParameter

    def test_root_length(self):
        """ tests if numerical root length agrees with analytic solutions at 4 points in time with two scales of dt"""
        self.root_example_rrp()
        times = np.array([0., 7., 15., 30., 60.])
        dt = np.diff(times)
        k = self.root.param().getK()  # maximal root length
        self.assertAlmostEqual(k, 100, 12, "example root has wrong maximal length")
        l = rootLength(times[1:], self.p0.r, k)  # analytical root length
        root = self.root.copy(self.plant)
        self.root_length_test(dt, l, 1)  # large dt
        self.root = root
        self.root_length_test(dt, l, 1000)  # very fine dt

    def test_root_length_including_laterals(self):
        """ tests if numerical root length agrees with analytic solution including laterals """
        self.root_example_rrp()
        times = np.array([0., 7., 15., 30., 60.])
        dt = np.diff(times)
        rp = self.root.getRootRandomParameter()
        p = self.root.param()  # rename
        k = p.getK()
        et = np.zeros((p.nob()))
        l = 0
        et[0] = rootAge(p.la - rp.ln / 2. + p.lb + l, p.r, k)
        for i in range(0, p.nob() - 1):  # calculate lateral emergence times
            l += p.ln[i]
            et[i + 1] = rootAge(p.la - rp.ln / 2. + p.lb + l, p.r, k)
        l = rootLength(times[1:], p.r, k)  # zero order lengths
        l1 = []
        r2 = self.p1.r
        k2 = self.p1.lmax  # consists of lateral zone only

        for t in times[1:]:
            l1.append(rootLateralLength(t, et, r2, k2))

        analytic_total = l + l1

        for subDX in [1, 1000]:
            numeric_total = []
            for t in times[1:]:
                root = self.root.copy(self.plant)
                self.root_length_test([t], [rootLength(t, p.r, k)], subDX)
                organs = self.root.getOrgans()
                nl = 0
                for o in organs:
                    nl += o.getParameter("length")
                numeric_total.append(nl);
                self.root = root
            for i in range(0, len(times[1:])):
                self.assertAlmostEqual(numeric_total[i], analytic_total[i], 10, "numeric and analytic total lengths do not agree in time step " + str(i + 1))

    def test_parameter(self):
        """ tests some parameters on sequential organ list """
        self.root_example_rrp()
        simtime = 30.
        self.root.simulate(simtime, False)
        organs = self.root.getOrgans()
        type, age, radius, order, ct = [], [], [], [], []
        for o in organs:
            type.append(o.getParameter("subType"))
            age.append(o.getParameter("age"))
            ct.append(o.getParameter("creationTime"))
            radius.append(o.getParameter("radius"))
            order.append(o.getParameter("order"))

        nol = round(self.root.getParameter("numberOfLaterals"))
        type_ = [1.]
        type_.extend([2.] * nol)
        self.assertEqual(type, type_, "getParameter: unexpected root sub types")
        self.assertEqual(order, type_, "getParameter: unexpected root order")  # +1, because of artificial parent root
        for i in range(0, nol):
            # print(i, nol, age[i], simtime - ct[i])
            self.assertAlmostEqual(age[i], simtime - ct[i], 10, "getParameter: age and creation time does not agree")

    def test_dynamics(self):
        """ tests if nodes created in last time step are correct """
        self.root_example_rrp()
        r = self.root
        r.simulate(.5, False)
        self.assertEqual(r.hasMoved(), False, "dynamics: node is creaetd during first step")
        r.simulate(1e-1, False)
        self.assertEqual(r.hasMoved(), True, "dynamics: node was expected to move, but did not")
        non = r.getNumberOfNodes()
        r.simulate(2.4, False)
        self.assertEqual(r.getOldNumberOfNodes(), non, "dynamics: wrong number of old nodes")
        dx = r.getRootRandomParameter().dx
        self.assertEqual(r.getNumberOfNodes() - non, round(2.4 * r.param().r / dx), "dynamics: unexpected number of new nodes")  # initially, close to linear growth

    def root_example_rrp2(self):
        """ an example used in the tests below, a main root with laterals """
        self.plant = pb.RootSystem()  # store organism (not owned by Organ, or OrganRandomParameter)
        p0 = pb.RootRandomParameter(self.plant)
        p0.name, p0.subType, p0.la, p0.lb, p0.lmax, p0.ln, p0.lnk, p0.r, p0.dx, p0.dxMin = "taproot", 1, 0.95, 0.8, 10., 1.05, 0.01, 0.8, 0.25, 0.2
        p0.successor = [[2]]
        p0.successorP = [[1.]]
        p1 = pb.RootRandomParameter(self.plant)
        p1.name, p1.subType, p1.lmax, p1.r, p1.dx = "lateral", 2, 2., 2., 2.

        self.plant.setOrganRandomParameter(p0)  # the organism manages the type parameters and takes ownership
        self.plant.setOrganRandomParameter(p1)
        srp = pb.SeedRandomParameter(self.plant)
        self.plant.setOrganRandomParameter(srp)

        print("root p0, initial parameters: lmax = ", p0.lmax, ", lb = ", p0.lb, ", la = ", p0.la, ", ln = ", p0.ln)
        param0 = p0.realize()  # set up root by hand (without a root system)
        print("root p0, realized parameters: lmax = ", sum((sum(param0.ln), param0.lb, param0.la)), ", lb = ", param0.lb, ", la = ", param0.la, ", mean ln = ", np.mean(param0.ln))
        if((param0.lb % p0.dx > 0) and (param0.lb % p0.dx < p0.dxMin * 0.99)):
            print("lb value does not fit with dx and dxMin")
            print(param0.lb % p0.dx)
        if((param0.la % p0.dx > 0) and (param0.la % p0.dx < p0.dxMin * 0.99)):
            print("la value does not fit with dx and dxMin")
            print(param0.la % p0.dx)
        if(any([(lni % p0.dx > 0 and  lni % p0.dx < p0.dxMin * 0.99) for lni in param0.ln])):
            print("ln value does not fit with dx and dxMin")

        param0.la, param0.lb = 0, 0  # its important parent has zero length, otherwise creation times are messed up
        parentroot = pb.Root(1, param0, True, True, 0., 0., pb.Vector3d(0, 0, -1), 0, 0, False, 0)  # takes ownership of param0
        parentroot.setOrganism(self.plant)
        parentroot.addNode(pb.Vector3d(0, 0, -1), 0)  # there is no nullptr in Python

        self.parentroot = parentroot  # store parent (not owned by child Organ)
        self.root = pb.Root(self.plant, p0.subType, pb.Vector3d(0, 0, -1), 0, self.parentroot , 0, 0)
        self.root.setOrganism(self.plant)
        self.p0 = p0

    def root_dxMin_test(self, dt):
        """ simulates a single root and checks length against analytic length """
        self.root_example_rrp2()
        nl, nl_th, non = [], [], []
        tot_dt = 0
        k = self.root.param().getK()  # maximal root length
        lb = self.root.param().lb
        la = self.root.param().la
        effectiveLa = la - np.mean(self.root.param().ln) / 2
        ln = np.concatenate((np.array([lb]), np.array(self.root.param().ln)))
        ln = np.cumsum(ln)
        for t in dt:
            self.root.simulate(t , True)
            tot_dt += t
            nl.append(self.root.getParameter("length"))
            l_th = rootLength(tot_dt, self.p0.r, k)  # analytical root length
            res = l_th % self.p0.dx

            if(res < self.p0.dxMin * 0.99):
                l_th -= res
            nl_th.append(l_th)
            non.append(self.root.getNumberOfNodes())
            # length from geometry
            poly = np.zeros((non[-1], 3))  #
            for i in range(0, non[-1]):
                v = self.root.getNode(i)
                poly[i, 0] = v.x
                poly[i, 1] = v.y
                poly[i, 2] = v.z
            d = np.diff(poly, axis = 0)
            length_segments = np.array([round(norm(di), 6) for di in d])
            if(np.min(length_segments) < self.p0.dxMin * 0.99):
                print("minimum segment length ", np.min(length_segments), " below dxMin ", self.p0.dxMin)
            if(self.p0.dx < np.max(length_segments) * 0.99):
                print("maximum segment length ", np.max(length_segments), " above dx ", self.p0.dx)
            if(len(np.where(ln <= (l_th - effectiveLa))[0]) != self.root.getParameter("numberOfLaterals")):
                print("\n\n\nnumeric and analytic number of laterals different: ", ln, l_th , effectiveLa, l_th - effectiveLa, len(np.where(ln <= (l_th - effectiveLa))[0]), self.root.getParameter("numberOfLaterals"))

        for i in range(0, len(dt)):
            self.assertAlmostEqual(nl_th[i], nl[i], 10, "numeric and analytic lengths do not agree in time step " + str(i + 1))
        print("end of test")

    def root_example_rtp2(self,delay_definition = 1):
        """ an example used in the tests below, a main root with laterals """
        self.partialiheading = pb.Vector3d.rotAB(0, 0)
        self.plant = pb.RootSystem()  # store organism (not owned by Organ, or OrganRandomParameter)
        p0 = pb.RootRandomParameter(self.plant)
        self.lmax_th = 100
        p0.name, p0.subType, p0.la, p0.lb, p0.lmax, p0.ln, p0.r, p0.dx, p0.dxMin = "main", 1, 10., 10., self.lmax_th, 1., 1.5, 1, 0.5
        p0.ldelay = 1.
        p0.successor = [[5]]
        p0.successorP = [[1.]]

        p1 = pb.RootRandomParameter(self.plant)
        p1.name, p1.subType, p1.lmax, p1.r, p1.dx, p1.dxMin = "lateral", 5, 5., 2., 1, 0.5
        p1.ldelay = 2.
        self.p0, self.p1 = p0, p1  # needed at later point
        self.plant.setOrganRandomParameter(p0)  # the organism manages the type parameters and takes ownership
        self.plant.setOrganRandomParameter(p1)

        srp = pb.SeedRandomParameter(self.plant)
        srp.delayDefinition = delay_definition #organ carries the delay of its laterals
        self.plant.setOrganRandomParameter(srp)
        # test == True => no need to give root parameter
        self.plant.initialize(verbose = False)
        paramS = srp.realize()
        self.seed = self.plant.getSeed()  #

        param0 = p0.realize()  # set up root by hand (without a root system)
        param0.la, param0.lb = 0, 0  # its important parent has zero length, otherwise creation times are messed up
        self.ons =pb.Vector3d(0., 0., 1.)
        parentroot = pb.Root(1, param0, True, True, 0., 0., self.partialiheading, 0, False, 0)  # takes ownership of param0
        parentroot.setOrganism(self.plant)
        parentroot.setParent(self.seed)
        parentroot.addNode(pb.Vector3d(0, 0, -3), 0)  # there is no nullptr in Python
        self.parentroot = parentroot  # store parent (not owned by child Organ)
        self.seed.addChild(self.parentroot)
        self.root = pb.Root(self.plant, p0.subType, 0, self.parentroot , 0)
        self.parentroot.addChild(self.root)
        self.root.setOrganism(self.plant)


    def test_new_delay_types(self):
        self.root_example_rtp2(delay_definition = 0) #depends on la
        r = self.root
        time = 100
        r.simulate(time, False)
        meanLn = self.root.getParameter("lnMean")# mean inter-lateral distance
        effectiveLa = max(self.root.getParameter("la")-meanLn/2, 0.)#effective apical distance, observed apical distance is in la-ln/2, la+ln/2
        
        for i in range(self.root.getNumberOfChildren()):  
            r = self.root.getChild(i)   
            rparent = r.getParent()
            lengthAtCreation = rparent.getLength(r.parentNI)
            ageLN = rparent.calcAge(lengthAtCreation)# age of root when lateral node is created
            ageLG = rparent.calcAge(lengthAtCreation+effectiveLa)# age of the root, when the lateral starts growing (i.e when the apical zone is developed)
            dl = ageLG-ageLN;   
            et = ageLN + dl
            self.assertAlmostEqual(r.getAge(), (time - et), 10, "numeric and analytic age of root n#" + str(i + 1) + " do not agree")
                
        self.root_example_rtp2(delay_definition = 1) #depends on ldelay of parent organ
        r = self.root
        time = 10
        r.simulate(time, False)
        for i in range(self.root.getNumberOfChildren()):  
            r = self.root.getChild(i)
            rsp = r.getParam()
            dl = r.getParent().getOrganRandomParameter().ldelay # only works because deviation == 0
            et = r.getParent().getNodeCT(r.parentNI) + dl
            self.assertAlmostEqual(r.getAge(), (time - et), 10, "numeric and analytic age of root n#" + str(i + 1) + " do not agree")
            
        self.root_example_rtp2(delay_definition = 2) #depends on ldelay of lateral organ
        r = self.root
        time = 10
        r.simulate(time, False)
        for i in range(self.root.getNumberOfChildren()):  
            r = self.root.getChild(i)
            rsp = r.getParam()
            dl = r.getOrganRandomParameter().ldelay # only works because deviation == 0
            et = r.getParent().getNodeCT(r.parentNI) + dl
            self.assertAlmostEqual(r.getAge(), (time - et), 10, "numeric and analytic age of root n#" + str(i + 1) + " do not agree")
if __name__ == '__main__':
    unittest.main()
