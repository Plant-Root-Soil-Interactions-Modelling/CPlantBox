import sys; sys.path.append("../.."); sys.path.append("../../src/")

import rsml.rsml_reader as rsml

import numpy as np
from scipy.optimize import minimize


class Root:
    """ A structure for easy analysis considering multiple measurements, 
        use parse_rsml to obtain a dictionary of Root objects with root ids as keys, 
        use initialize_roots to create a linked list (parents and laterals) """

    def __init__(self):
        """ creates an empty Root """
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
        self.r = -1.  # initial growth rate (of root type)
        self.k = -1.  # maximal root length (of root type)
        self.ldelay = -1.  # in case of delay based growth (use RootSystem::intitializeDB())

    def add_measurement(self, time:float, i:int, poly:np.array, prop:dict, fun:dict, roots:dict):
        """ creates the i-th root from the polylines and adds it to the roots dictionary """
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

    def initialize(self, roots:dict):
        """ connects the roots creating parent and laterals expects prop['parent-poly'], prop['parent-node'] """
        t0 = self.measurement_times[0]
        parent_id = self.props[t0]('parent-poly')
        if parent_id > -1:
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
        """ sorts the laterals after their parent node index, call initialize for all roots first """
        for t in self.measurement_times:
            pni = np.array([l.parent_node for l in self.laterals.setdefault(t, [])])
            ii = np.argsort(pni)
            l_ = []
            for i in ii:
                 l_.append(self.laterals[t][i])
            self.laterals[t] = l_

    def length(self, i0:int = 0, iend:int = -1, t:float = -1.):
        """ length between two indices @param i0 and @param iend at measurement @param m """
        if t < 0:
            t = self.measurement_times[-1]
        if iend == -1:
            iend = self.nodes[t].shape[0] - 1
        l = 0.
        if i0 < iend:
            for i in range(int(i0), int(iend)):
                l += np.linalg.norm(self.nodes[t][i] - self.nodes[t][i + 1])
        return l

    def order(self, i:int = 0):
        """ recursively determines root order """
        assert self.id > -1, "Root.order: uninitialized root"
        if self.parent:
            return self.parent.order(i + 1)
        else:
            return i

    def set_emergence_time(self, et:float):
        """ sets the emergance times, calculates the root ages """
        self.emergence_time = et
        for t in self.measurement_times:
            self.ages[t] = max(t - et, 0.)

    def calc_growth_rate(self, r:float, k:float):
        """ sets the maximal root length [cm] and computes the individual root initial growth rate self.r [cm/day], 
        and lateral delay time (based on the apical length).  The caluclation is based on self.emergence_time, i.e. call set_emergence_time first
        Furthermore sets the emergence times of its laterals. """
        self.k = k
        r_, ldelay_ = [], []
        for t in self.measurement_times:
            l = self.length(0, -1, t)
            age = self.ages[t]
            if l > 0 and age > 0:
                r_.append(-k / age * np.log(1 - l / k))
            else:
                r_.append(np.nan)
        self.r = np.nanmean(np.array(r_))
        if np.isnan(self.r):
            rr = r
        else:
            rr = self.r  # if we have something we can use...
        for t in self.measurement_times:
            if self.laterals[t]:
                length_la = self.length(0, self.laterals[t][-1].parent_node, t)  # length - la
                ld = self.ages[t] - Root.negexp_age(length_la + self.la, rr, self.k)
                ldelay_.append(max(ld, 0.))
            else:
                ldelay_.append(np.nan)
        self.ldelay = np.nanmean(np.array(ldelay_))
        lm = self.measurement_times[-1]  # last measurement
        # print("\nRoot ", self.id)
        for i, l in enumerate(self.laterals[lm]):
            tl = self.length(0, l.parent_node, -1) + self.la
            lateral_et = Root.negexp_age(tl, rr, self.k) + self.emergence_time
            l.set_emergence_time(lateral_et)
            # print(lateral_et, l.ages[lm])
        # print("")

    def calc_params(self):
        """ retrieves la, lb, ln, theta, a """
        lm = self.measurement_times[-1]  # last measurement
        la_, lb_, ln_, a_ = [], [], [], []
        for t in self.measurement_times:
            l = self.laterals[t]
            if l:
                la_.append(self.length(l[-1].parent_node, -1, t))
                lb_.append(self.length(0, l[0].parent_node, t))
            else:
                la_.append(np.nan)
                lb_.append(np.nan)
        self.la = np.mean(np.array(la_))
        self.lb = np.nanmean(np.array(lb_))
        l = self.laterals[lm]  # ln is calculated from the last measurement
        if l:
            for i in range(0, len(l) - 1):
                ln = self.length(l[i].parent_node, l[i + 1].parent_node, t)
                ln_.append(ln)
            if len(l) - 1 == 0:
                ln_.append(np.nan)
        else:
            ln_.append(np.nan)
        self.ln = np.nanmean(np.array(ln_))
        if np.isnan(self.ln):
            for i in range(0, len(l) - 1):
                print(self.length(l[i].parent_node, l[i + 1].parent_node, t))
        # print(self.la, self.lb, self.ln)
        if self.nodes[lm].shape[0] > 1:
            v1 = self.nodes[lm][1] - self.nodes[lm][0]  # theta is calculated from the last measurement
            v1 = v1 / np.linalg.norm(v1)
            v2 = np.array([0., 0., -1])
            if self.parent:
                pni = int(self.parent_node)
                if pni > 1:
                    v2 = self.parent.nodes[lm][pni] - self.parent.nodes[lm][pni - 2]  # Kutschera file seems to store nodes twice
                    v2 = v2 / np.linalg.norm(v2)
            self.theta = np.arccos(np.clip(np.dot(v1, v2), -1.0, 1.0))
        else:
            self.theta = np.nan
        # print(self.theta)
        for t in self.measurement_times:
            n = self.nodes[t].shape[0]
            a = 0.
            for i in range(0, n):
                a += self.funs[t]('diameter', i) / n / 2.
            a_.append(a)
        self.a = np.mean(np.array(a_))

    def negexp_length(t, r, k):
        """ returns the root length [cm] assuming negative exponential growth 
        @param t root age [day] 
        @param r initial root growth [cm/day] 
        @param k maximal root length [cm] """
        return k * (1 - np.exp(-(r / k) * t))

    def negexp_age(l, r, k):
        """ returns the root age [day] assuming negative exponential growth, 
        @param l root length [cm] 
        @param r initial root growth [cm/day] 
        @param k maximal root length [cm] """
        return -k / r * np.log(1 - l / k)


def connect(node, base_polyline):
    """ connects the node to the closest point in base_polyline """
    min_dist = 1.e8  # much
    for i, n in enumerate(base_polyline):
        dist = np.linalg.norm(n - node)
        if dist < min_dist:
            min_dist = dist
            min_i = i
    return min_i


def parse_rsml(rsml_name:list, time:list):
    """ creates a root dictionary from a rsml file
    @param rsml_name     file name 
    @param time          measurement time  """
    roots = {}
    poly, prop, fun, _ = rsml.read_rsml(rsml_name)
    assert len(poly) == len(prop['parent-poly']), "parse_rsml: wrong number of parent-poly tags"
    if not 'parent-node' in prop:  # reconstruct...
        print("parse_rsml: no parent-node tag found, reconstructing...")
        prop['parent-node'] = []
        for i in range(0, len(poly)):
            if prop['parent-poly'][i] >= 0:
                pni = connect(np.array(poly[i][0]), poly[prop['parent-poly'][i]])
            else:
                pni = -1
            # print("parent", prop['parent-poly'][i], "pni", pni)
            prop['parent-node'].append(pni)
    nopn = len(prop['parent-node'])
    assert len(poly) == nopn, "parse_rsml: wrong number of parent-node tags"
    for i in range(0, len(poly)):
        r = Root()
        r.add_measurement(time, i, poly, prop, fun, roots)
    initialize_roots(roots)
    return roots


def parse_rsmls(rsml_names:list, times:list):
    """ creates a root dictionary from a rsml file
    @param rsml_names    file names 
    @param times         measurement times  """
    assert len(rsml_names) == len(times), "parse_rsmls: number of file names must equal the number of measurement times"
    root_list = []
    for i, n in enumerate(rsml_names):
        root_list.append(parse_rsml(n, times[i]))
    return merge_measurements(root_list)


def initialize_roots(roots:dict):
    """ calls initialize for each root, and sorts the laterals by their parent node index, 
    @param roots: root dicitionaray """
    for r in roots.values():
        r.initialize(roots)
    for r in roots.values():
        r.sort_laterals()


def get_order(i:int, roots:dict):
    """ obtain a list of order @param i roots from the root dicitionaray @param roots """
    r_ = []
    for r in roots.items():
        if r[1]:
            if r[1].order() == i:
                r_.append(r[1])
    return r_


def get_params(roots:list, time:float):
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
    return np.array([[np.nanmean(p), np.nanstd(p)] for p in params])


def merge_measurements(root_list:list):
    """ merges multiplie measurements of the same plant type, the roots must be initialized """
    roots = {}
    for i, r in enumerate(root_list):
        for r_ in r.values():
            r_.id = r_.id + i * 10000
            roots[r_.id] = r_
    return roots


#
# Estimation of parameters
#
def get_see_rk(lengths:np.array, ages:np.array, r:float, k:float):
    """ returns the strandard error of the estimate for r and k (under questionable conditions) """
    assert lengths.shape == ages.shape, "the number of measured lengths must equal measured ages"
    r_ = lambda l, age:-k / age * np.log(1 - l / k)
    k_ = lambda l, age: l / (1 - np.exp(-r / k * age))  # still dependent on k, but use as approximation for now
    mr_ = []  # r by measured values l and age for a fixed k
    mk_ = []  # k by measured values l and age for a fixed r
    for i, l in enumerate(lengths):
        age = ages[i]
        if l < k and age > 0:
            mr_.append(r_(l, age))
        if age > 0:
            mk_.append(k_(l, age))
    mr_ = np.array(mr_)
    mk_ = np.array(mk_)
    see_r = np.sqrt(np.sum(np.square(mr_ - r * np.ones(mr_.shape))) / (mr_.shape[0] - 2))
    see_k = np.sqrt(np.sum(np.square(mk_ - k * np.ones(mk_.shape))) / (mk_.shape[0] - 2))
    return see_r, see_k


def estimate_set_order0_rate(roots:dict, lmax:float, time:float):
    """ estimates """
    basals = get_order(0, roots)
    basal_ids = np.array([r.id for r in basals], dtype = np.int64)
    basal_lengths = np.array([r.length() for r in basals])
    basal_ages = np.array([time] * len(basals))
    ii = np.argsort(basal_lengths)  # sort by ascending lengths
    basal_ids = basal_ids[ii]
    basal_lengths = basal_lengths[ii]
    basal_ages = basal_ages[ii]
    res, f, ages = estimate_order0_rrate(basal_lengths, lmax, time, 1.5)
    rate, r = res.x[0], res.x[1]
    rs, _ = get_see_rk(basal_lengths, ages, r, lmax)
    ages = np.zeros(basal_ids.shape)
    for i, _ in enumerate(basal_ids):
        ages[basal_ids.shape[0] - i - 1] = max(time - i * rate, 0.)
    for i, id in enumerate(basal_ids):
        et = time - ages[i]
        roots[id].set_emergence_time(et)
    for id in basal_ids:  # set order 0 growth rate and maximal length
        roots[id].calc_growth_rate(r, lmax)
    return rate, r, rs


def estiamte_emergance_order0(lengths:np.array, ages:np.array, r:float, k:float):
    """ fits the emergance time of 2nd basal roots """
    f = lambda x: target_length(r, k, lengths, ages - np.ones(ages.shape) * x[0])
    x0 = [np.max(ages)]
    res = minimize(f, x0, method = 'Nelder-Mead', tol = 1e-6)
    return res, f


def estimate_order0_rate(lengths:np.array, r:float, k:float, time:float):
    """ fits basal prodcution rate [day-1] for given initial growth rate and maximal root length
    @param lengths list of root lengths [cm] 
    @param r initial root length [cm] 
    @param k maximal root length [cm]
    @param time maximal measurement time """
    f = lambda x: target_rate(x[0], lengths, r, k, time)
    x0 = [time / lengths.shape[0]]
    res = minimize(f, x0, method = 'Nelder-Mead', tol = 1e-6)  # bounds and constraints are possible, but method dependent
    n = lengths.shape[0]
    ages = np.zeros(lengths.shape)
    for i in range(0, n):
        ages[n - i - 1] = max(time - i * res.x[0], 0.)
    return res, f, ages  # e.g. res.x[0], f(res.x[0])


def estimate_order0_rrate(lengths:np.array, k:float, time:float, r0:float = 1.):
    """ fits basal prodcution rate [day-1] anjnd initial growth rate and maximal root length, 
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


def estimate_r(lengths:np.array, ages:np.array, k:float):
    """ fits initial growth rate r [cm/day], assumes maximal root lenght k as fixed (e.g. by literature value) """
    assert len(lengths) == len(ages), "estimate_r: number of root lengths must equal number of root ages"
    f = lambda x0: target_length(x0, k, lengths, ages)
    x0 = [1.]  # initial value
    res = minimize(f, x0, method = 'Nelder-Mead', tol = 1e-6)  # bounds and constraints are possible, but method dependent
    return res, f  # e.g. res.x[0], f(res.x[0])


def estimate_rk(lengths:np.array, ages:np.array):
    """ fits initial growth rate r, and maximal root lenght k """
    assert lengths.shape == ages.shape, "estimate_rk: number of root lengths must equal number of root ages"
    f = lambda x0: target_length(x0[0], x0[1], lengths, ages)
    x0 = [5., 200]  # initial value
    res = minimize(f, x0, method = 'Nelder-Mead', tol = 1e-6)  # bounds and constraints are possible, but method dependent
    return res, f  # e.g. res.x[0], res.x[1], f(res.x)


def target_length(r:float, k:float, lengths:np.array, ages:np.array):
    """ target function for optimization root target length [cm],
    @param r initial root length [cm], @param k maximal root length [cm]
    @param lengths root lengths [cm], @param ages root ages [day] corresponding to root lengths"""
    assert lengths.shape == ages.shape, "target_length: number of root lengths must equal number of root ages"
    sum = 0.
    for i in range(0, lengths.shape[0]):
        l = Root.negexp_length(ages[i], r, k)
        sum += (l - lengths[i]) ** 2
    return np.sqrt(sum / (lengths.shape[0] - 2))  # -(k+1)


def target_length2(r:float, lmax:float, lengths:np.array, ages:np.array, time:float):
    """ target function for optimization root target length [cm] minimizing the distance between curve and point. 
    @param r initial root length [cm], @param k maximal root length [cm]
    @param lengths root lengths [cm], @param ages root ages [day] corresponding to root lengths"""
    assert lengths.shape == ages.shape, "target_length: number of root lengths must equal number of root ages"
    sum = 0.
    for i in range(0, lengths.shape[0]):
        f = lambda t_: (Root.negexp_length(t_, r, lmax) - lengths[i]) ** 2 + (t_ - ages[i]) ** 2  # distance between curve and (lengths[i], ages[i])
        res = minimize(f, [ages[i]], method = 'Nelder-Mead', tol = 1e-6)
        t_ = res.x[0]
        sum += (Root.negexp_length(t_, r, lmax) - lengths[i]) ** 2 + (t_ - ages[i]) ** 2
    return  np.sqrt(sum / (lengths.shape[0] - 2))


def target_rate(rate:float, lengths:np.array, r:float, lmax:float, time:float):
    """ target function for estimating the linear base root production rate [day-1],
    @param rate: linear base root production rate (delay between basal root emergence) [day-1] 
    @param lengths: basal root lengths as numpy array sorted ascending [cm]
    @param r: initial growth rate [cm/day]
    @param lmax: maximal root length [cm] 
    @param time: maximal measurement time
    @return error 
     """
    rate = max(rate, 0.)
    n = lengths.shape[0]
    ages = np.zeros(lengths.shape)
    for i in range(0, n):
        ages[n - i - 1] = max(time - i * rate, 0.)
    x = target_length(r, lmax, lengths, ages)
    # x = target_length2(r, k, lengths, ages, time)
    return x

