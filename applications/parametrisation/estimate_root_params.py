import sys;  sys.path.append("../..")
import numpy as np
from scipy.optimize import minimize

import rsml_reader as rsml

import plantbox as pb


class Root:
    """ A structure for easy analysis considering multiple measurements, 
        use parse_rsml to obtain a dictionary of Root objects with root ids as keys, 
        use initialize_roots to create linked list (parents and laterals) """

    def __init__(self):
        """ creates empty Root """
        self.id = -1
        self.measurements = 0  # number of measurements
        self.measurement_times = []
        # parent
        self.parent = None  # self.parent.id
        self.parent_node = -1
        self.parent_base_length = -1
        # laterals
        self.laterals = {}
        self.nodes = {}
        # RSML additional data
        self.props = {}
        self.funs = {}
        # parameters
        self.la = 0.  # of this root
        self.lb = 0.  # of this root
        self.ln = 0.  # of this root
        self.a = 0.  # mean radius of this root
        self.theta = 0.  # of this root
        # must be set from outside
        self.emergence_time = 0.  #
        self.ages = {}
        self.r = 0.  # initial growth rate (of root type)
        self.k = 0.  # maximal root length (of root type)

    def add_measurement(self, time :float, i :int, poly :np.array, prop :dict, fun :dict, roots :dict):
        """ creates the i-th root from the polylines and adds it to the roots dictionary"""
        if self.measurements == 0:
            self.id = i
        assert self.id == i, "Root.add_measurement: root id changed"
        self.measurements += 1
        self.measurement_times.append(time)
        self.nodes[time] = np.array([np.array(n) for n in poly[i]])  # convert to 2d numpy array
        self.props[time] = (lambda key: prop[key][i])  # restrict to root i
        self.funs[time] = (lambda key, j: fun[key][i][j])  # restrict to root i
        assert self.id > -1, "Root.add_measurement: self.id has negative index"
        roots[self.id] = self

    def initialize(self, roots :dict):
        """ connects the roots creating parent and laterals expects prop['parent-poly'],  prop['parent-node']  """
        t0 = self.measurement_times[0]
        parent_id = self.props[t0]('parent-poly')
        if parent_id > 0:
            self.parent = roots[parent_id]
        else:
            self.parent = None
        self.parent_node = self.props[t0]('parent-node')
        if self.parent:
            self.parent_base_length = self.parent.length(0, self.parent_node)
            for t in self.measurement_times:
                self.parent.laterals.setdefault(t, []).append(self)  # register laterals
        else:
            self.parent_base_length = 0.
        for t in self.measurement_times[1:]:  # check integrity for following time steps
            if self.parent:
                assert self.parent.id == self.props[t]('parent-poly'), "Root.initialize: parent root id changed in measurements"
                assert self.parent_node == self.props[t]('parent-node'), "Root.initialize: parent node changed in measurements"

    def sort_laterals(self):
        """ sorts the laterals after their parent node index, call initialize() for all roots first """
        for t in self.measurement_times:
            l_ = []
            for i, l in enumerate(self.laterals.setdefault(t, [])):
                l_.append((l.parent_node, i))
            l_.sort()
            self.laterals[t] = [self.laterals[t][l_[j][1]] for j in range(0, len(l_))]

    def length(self, i0 :int = 0, iend :int = -1, t :float = -1.):
        """ length between two indices @param i0 and @param iend at measurement @param m """
        if t < 0:
            t = self.measurement_times[-1]
        if iend == -1:
            iend = self.nodes[t].shape[0]
        l = 0.
        for i in range(int(i0), int(iend) - 1):
            l += np.linalg.norm(self.nodes[t][i] - self.nodes[t][i + 1])
        return l

    def order(self, i :int = 0):
        """ recursively determines root order """
        assert self.id > -1, "Root.order: uninitialized root"
        if self.parent:
            return self.parent.order(i + 1)
        else:
            return i

    def set_emergence_time(self, et :float):
        """ sets the emergance times, calculates the root ages """
        self.emergence_time = et
        for t in self.measurement_times:
            self.ages[t] = max(t - et, 0.)

    def calc_growth_rate(self, r :float, k :float):
        """ sets the maximal root length [cm] and computes the root initial growth rate [cm/day], 
        call set_emergence_time first ! """
        self.k = k
        r_ = []
        for t in self.measurement_times:
            l = self.length(0, -1, t)
            age = self.ages[t]
            if l > 0 and age > 0:
                r_.append(-k / age * np.log(1 - l / k))
            else:
                r_.append(r)
        self.r = np.mean(np.array(r_))
        lm = self.measurement_times[-1]  # last measurement
        for i, l in enumerate(self.laterals[lm]):
            tl = self.length(0, l.parent_node, -1) + self.la
            lateral_age = Root.negexp_age(tl, self.r, self.k) + self.emergence_time
            l.set_emergence_time(lateral_age)

    def calc_params(self):
        """ retrieves la, lb, ln, theta, a """
        lm = self.measurement_times[-1]  # last measurement
        la_, lb_, ln_, a_ = [], [], [], []
        for t in self.measurement_times:
            l = self.laterals[t]
            if l:
                la_.append(self.length(0, l[0].parent_node, t))
                lb_.append(self.length(l[-1].parent_node, -1, t))
            else:
                la_.append(self.length(0, -1, t))
                lb_.append(0.)
        self.la = np.mean(np.array(la_))
        self.lb = np.mean(np.array(lb_))
        l = self.laterals[lm]  # ln is calculated from the last measurement
        if l:
            for i in range(0, len(l) - 1):
                ln_.append(self.length(l[i].parent_node, l[i + 1].parent_node, t))
            if len(l) - 1 == 0:
                ln_.append(0.)
        else :
            ln_.append(0.)
        self.ln = np.mean(np.array(ln_))
        if np.isnan(self.ln):
            print(np.array(ln_))
            for i in range(0, len(l) - 1):
                print(self.length(l[i].parent_node, l[i + 1].parent_node, t))
        # print(self.la, self.lb, self.ln)
        v1 = self.nodes[lm][1] - self.nodes[lm][0]  # theta is calculated from the last measurement
        v1 = v1 / np.linalg.norm(v1)
        v2 = np.array([0., 0., -1])
        if self.parent:
            pni = int(self.parent_node)
            if pni > 1:
                v2 = self.parent.nodes[lm][pni] - self.parent.nodes[lm][pni - 2]  # Kutschera file seems to store nodes twice
                v2 = v2 / np.linalg.norm(v2)
        self.theta = np.arccos(np.clip(np.dot(v1, v2), -1.0, 1.0))
        # print(self.theta)
        for t in self.measurement_times:
            n = self.nodes[t].shape[0]
            a = 0.
            for i in range(0, n):
                a += self.funs[t]('diameter', i) / n / 2.
            a_.append(a)
        self.a = np.mean(np.array(a_))

    def negexp_length(t, r, k):
        """ root length [cm] according to negative exponential growth, 
        @param t root age [day], @param r initial root growth [cm/day], @param k maximal root length [cm]"""
        return k * (1 - np.exp(-(r / k) * t))

    def negexp_age(l, r, k):
        """ root age [day] according to negative exponential growth, 
        @param l root length [cm], @param r initial root growth [cm/day], @param k maximal root length [cm]"""
        return -k / r * np.log(1 - l / k)


def initialize_roots(roots :dict):
    """ calls initialize for each root, and sorts the laterals by their parent node index, 
    @param roots root dicitionaray """
    for r in roots.items():
        if r[1]:
            r[1].initialize(roots)
    for r in roots.items():
        if r[1]:
            r[1].sort_laterals()


def parse_rsml(rsml_names :list, times :list):
    """ creates a root dictionary 
    @param rsml_names list of rsml names of multiple measurements of the same plant, 
    @param times measurement times  """
    assert len(rsml_names) == len(times), "parse_rsml: length of file name list and time list must be equal"
    roots = {}
    for j, n in enumerate(rsml_names):
        poly, prop, fun = rsml.read_rsml(n)
        nor = len(poly)
        nopp = len(prop['parent-poly'])
        nopn = len(prop['parent-node'])
        assert nor == nopp, "parse_rsml: wrong number of parent-poly tags"
        assert nor == nopn, "parse_rsml: wrong number of parent-node tags"
        for i in range(0, len(poly)):
            r = Root()
            r.add_measurement(times[j], i, poly, prop, fun, roots)
    initialize_roots(roots)
    return roots


def merge_plants(plants :list):
    pass


def get_order(i :int, roots :dict):
    """ obtain a list of order @param i roots from the root dicitionaray @param roots """
    r_ = []
    for r in roots.items():
        if r[1]:
            if r[1].order() == i:
                r_.append(r[1])
    return r_


def get_params(roots :list, time :float):
    """ takes mean and std of la, lb, ln, a, theta of measurement tiem @param time """
    # calculate la, lb, ln, a, theta
    for r in roots:
        r.calc_params()
    la = np.array([r.la for r in roots])
    lb = np.array([r.lb for r in roots])
    ln = np.array([r.ln for r in roots])
    a = np.array([r.a for r in roots])
    theta = np.array([r.theta for r in roots])
    params = [la, lb, ln, a, theta]
    return np.array([[np.mean(p), np.std(p)] for p in params])


#
# Estimation of parameters
#
def estimate_order0_rate(lengths :np.array, r :float, k :float, time :float):
    """ fits basal prodcution rate [day-1] for given initial growth rate and maximal root length, 
    @param lengths list of root lengths [cm], 
    @param r initial root length [cm], @param k maximal root length [cm], @param time maximal measurement time"""
    f = lambda x: target_rate(x[0], lengths, r, k, time)
    x0 = [time / lengths.shape[0]]
    res = minimize(f, x0, method = 'Nelder-Mead', tol = 1e-6)  # bounds and constraints are possible, but method dependent
    n = lengths.shape[0]
    ages = np.zeros(lengths.shape)
    for i in range(0, n):
        ages[n - i - 1] = max(time - i * res.x[0], 0.)
    return res, f, ages  # e.g. res.x[0], f(res.x[0])


def estimate_order0_rrate(lengths:np.array, k :float, time :float, r0 :float = 1.):
    """ fits basal prodcution rate [day-1] for given initial growth rate and maximal root length, 
    @param lengths list of root lengths [cm], 
    @param r initial root length [cm], @param k maximal root length [cm], @param time maximal measurement time"""
    f = lambda x: target_rate(x[0], lengths, x[1], k, time)
    x0 = [time / lengths.shape[0], r0]  # initial value
    res = minimize(f, x0, method = 'Nelder-Mead', tol = 1e-6)  # bounds and constraints are possible, but method dependent
    n = lengths.shape[0]
    ages = np.zeros(lengths.shape)
    for i in range(0, n):
        ages[n - i - 1] = max(time - i * res.x[0], 0.)
    return res, f, ages  # e.g. res.x[0], res.x[1], f(res.x)


def estimate_r(lengths :np.array, ages :np.array, k :float):
    """ fits initial growth rate r [cm/day], assumes maximal root lenght k as fixed (e.g. by literature value) """
    assert len(lengths) == len(ages), "estimate_r: number of root lengths must equal number of root ages"
    f = lambda x0: target_length(x0, k, lengths, ages)
    x0 = [1.]  # initial value
    res = minimize(f, x0, method = 'Nelder-Mead', tol = 1e-6)  # bounds and constraints are possible, but method dependent
    return res, f  # e.g. res.x[0], f(res.x[0])


def estimate_rk(lengths :np.array, ages :np.array):
    """ fits initial growth rate r, and maximal root lenght k """
    assert lengths.shape == ages.shape, "estimate_rk: number of root lengths must equal number of root ages"
    f = lambda x0: target_length(x0[0], x0[1], lengths, ages)
    x0 = [5., 200]  # initial value
    res = minimize(f, x0, method = 'Nelder-Mead', tol = 1e-6)  # bounds and constraints are possible, but method dependent
    return res, f  # e.g. res.x[0], res.x[1], f(res.x)


def target_length(r :float, k :float, lengths :np.array, ages :np.array):
    """ target function for optimization root target length [cm],
    @param r initial root length [cm], @param k maximal root length [cm]
    @param lengths root lengths [cm], @param ages root ages [day] corresponding to root lengths"""
    assert lengths.shape == ages.shape, "target_length: number of root lengths must equal number of root ages"
    sum = 0.
    for i in range(0, lengths.shape[0]):
        l = Root.negexp_length(ages[i], r, k)
        sum += (l - lengths[i]) ** 2
    return np.sqrt(sum)


def target_length2(r :float, k :float, lengths :np.array, ages :np.array, time :float):
    """ target function for optimization root target length [cm] using distance between curve and point. 
    @param r initial root length [cm], @param k maximal root length [cm]
    @param lengths root lengths [cm], @param ages root ages [day] corresponding to root lengths"""
    assert lengths.shape == ages.shape, "target_length: number of root lengths must equal number of root ages"
    sum = 0.
    for i in range(0, lengths.shape[0]):
        f = lambda t_: (Root.negexp_length(t_, r, k) - lengths[i]) ** 2 + (t_ - ages[i]) ** 2  # distance between curve and (lengths[i], ages[i])
        res = minimize(f, [ages[i]], method = 'Nelder-Mead', tol = 1e-6)
        t_ = res.x[0]
        sum += (Root.negexp_length(t_, r, k) - lengths[i]) ** 2 + (t_ - ages[i]) ** 2
    return  np.sqrt(sum)


def target_rate(rate :float, lengths :np.array, r :float, k :float, time :float):
    """ target function for optimization a linear base root production rate [day-1],
    @param rate linear base root production rate [day-1], @param  lengths root lengths [cm]
    @param r initial root length [cm], @param k maximal root length [cm], @param time maximal measurement time
     """
    rate = max(rate, 0.)
    n = lengths.shape[0]
    ages = np.zeros(lengths.shape)
    for i in range(0, n):
        ages[n - i - 1] = max(time - i * rate, 0.)
    x = target_length(r, k, lengths, ages)
    # x = target_length2(r, k, lengths, ages, time)
    return x

