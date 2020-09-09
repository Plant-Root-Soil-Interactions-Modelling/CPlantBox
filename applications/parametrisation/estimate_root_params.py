import sys;  sys.path.append("../..")
import numpy as np
from scipy.optimize import minimize
import rsml_reader as rsml

import plantbox as pb


class Root:
    """ A structure for easy analysis considering multiple measurements """

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
        self.parent = roots[parent_id]  # for -1 is None
        self.parent_node = self.props[t0]('parent-node')
        if self.parent:
            self.parent_base_length = self.parent.length(0, int(self.parent_node))
            for t in self.measurement_times:
                self.parent.laterals.setdefault(t, []).append(self)  # register laterals
        else:
            self.parent_base_length = 0.
        for t in self.measurement_times[1:]:  # check integrity for following time steps
            if self.parent:
                assert self.parent.id == self.props[t]('parent-poly'), "Root.initialize: parent root id changed in measurements"
                assert self.parent_node == self.props[t]('parent-node'), "Root.initialize: parent node changed in measurements"

    def sort_laterals(self):
        """ sorts the laterals after their parent node index """
        for t in self.measurement_times:
            l_ = []
            for i, l in enumerate(self.laterals.setdefault(t, [])):
                l_.append((l.parent_node, i))
                l_.sort()
            self.laterals[t] = [self.laterals[t][l_[j][1]] for j in range(0, len(l_))]
#             pns = [l_[j][0] for j in range(0, len(l_))]
#             print(pns)
#             print([l.id for l in self.laterals[t]])

    def length(self, i0 :int = 0, iend :int = -1, t :float = -1.):
        """ length between two indices @param i0 and @param iend at measurement @param m """
        if t < 0:
            t = self.measurement_times[-1]
        if iend == -1:
            iend = self.nodes[t].shape[0]
        l = 0.
        for i in range(i0, iend - 1):
            l += np.linalg.norm(self.nodes[t][i] - self.nodes[t][i + 1])
        return l

    def order(self, i :int = 0):
        """ recursively determines root order """
        assert self.id > -1, "Root.order: uninitialized root"
        if self.parent:
            return self.parent.order(i + 1)
        else:
            return i


def initialize_roots(roots):
    """ call initialize for each root, and sort the laterals """
    roots[-1] = None
    for r in roots.items():
        if r[1]:
            r[1].initialize(roots)
    for r in roots.items():
        if r[1]:
            r[1].sort_laterals()


def parse_rsml(rsml_names, times):
    """ creates the roots dictionary """
    assert len(rsml_names) == len(times), "parse_rsml: length of file name list and time list must be equal"
    roots = {}
    for j, n in enumerate(rsml_names):
        poly, prop, fun = rsml.read_rsml(n)
        # use 'parent-poly', 'parent-node']
        nor = len(poly)
        nopp = len(prop['parent-poly'])
        nopn = len(prop['parent-node'])
        assert nor == nopp, "parse_rsml: wrong number of parent-poly tags"
        assert nor == nopn, "parse_rsml: wrong number of parent-node tags"
        for i in range(0, len(poly)):
            r = Root()
            r.add_measurement(times[j], i, poly, prop, fun, roots)
    initialize_roots(roots)


def negexp_length(t, r, k):
    return k * (1 - np.exp(-(r / k) * t))


def negexp_age(l, r, k):
    return -k / r * np.log(1 - l / k)


def target(r, k, length, times):
    """ target function for optimization for fit_taproot_r, fit_taproot_rk """
    sum = 0.
    for i, t in enumerate(times):
        sim_l = negexp_length(t, r, k)
        sum_t = 0.
        if isinstance(length[i], float):
            sum_t += (sim_l - length[i]) ** 2
        else:
            for l in length[i]:
                sum_t += (sim_l - l) ** 2

        sum += np.sqrt(sum_t)
    return sum


def fit_taproot_r(length, times, k):
    """ fits initial growth rate r, assumes maximal root lenght k as fixed (e.g. literature value) """
    assert(len(length) == len(times))
    f = lambda x0: target(x0, k, length, times)
    x0 = [1.]
    res = minimize(f, x0, method = 'Nelder-Mead', tol = 1e-6)  # bounds and constraints are possible, but method dependent
    return res.x[0], f(res.x[0])


def fit_taproot_rk(length, times):
    """ fits initial growth rate r, and maximal root lenght k """
    assert(len(length) == len(times))
    f = lambda x0: target(x0[0], x0[1], length, times)
    x0 = [5., 200]
    res = minimize(f, x0, method = 'Nelder-Mead', tol = 1e-6)  # bounds and constraints are possible, but method dependent
    return res.x[0], res.x[1], f(res.x)


def target2(delay, times, numbers, i_n):
    """ target function for optimization for fit_number_of_roots"""
    sum = 0.
    for i, t in enumerate(times):
        sim_n = i_n + t / delay  # missing round (continious seems easier to optimize)
        s = 0.
        for j in range(0, len(numbers[i])):
            s += (numbers[i][j] - sim_n) ** 2
        sum += np.sqrt(s)
    return sum


def fit_number_of_roots(times, numbers, initial_number):
    """ we fit the delay between emergence of seminals with, 
        number_of_roots = inital_number + round(t/delay)
    """
    assert(len(numbers) == len(times))
    f = lambda x0 : target2(x0, times, numbers, initial_number)
    x0 = [1.]  # days
    res = minimize(f, x0, method = 'Nelder-Mead', tol = 1e-6)  # linear regression would be enough in this case
    return res.x[0]

# def fit_seminal_rk(length, times):
#     """ fits initial growth rate r, and maximal root lenght k """
#     assert(len(length.shape[0] == len(times))
#     f = lambda x0: target(x0[0], x0[1], length, times)
#     x0 = [5., 200]
#     res = minimize(f, x0, method='Nelder-Mead', tol=1e-6)  # bounds and constraints are possible, but method dependent
#     # print(x)
#     return res.x[0], res.x[1]


def connect(node, base_polyline):
    """
    connects the node to the closest point in base_polyline
    """
    min_dist = 1.e8  # much
    for i, n in enumerate(base_polyline):
        dist = np.linalg.norm(n - node)
        if dist < min_dist:
            min_dist = dist
            min_i = i
    return min_i, min_dist


def reconstruct_laterals(polylines, properties, base_polyline, snap_radius = 0.5):
    """
    connect all lateal roots to the base root (all roots other than tap need will have order 1)
    """
    s = 0
    npl, nprop = [], {}
    pp = properties["parent-poly"]  # generated by read_rsml
    pn = properties["parent-node"]  # generated by read_rsml
    for i, pi in enumerate(pp):
        if pi == -1:  # tap root
            npl.append(polylines[i])
            for k in properties.keys():
                nprop.setdefault(k, []).append(properties[k][i])
        elif pi == 0:  # lateral root
            npl.append(polylines[i])
            for k in properties.keys():
                nprop.setdefault(k, []).append(properties[k][i])
        elif pi > 0:  # attach to base root
            n = np.array(polylines[i][0])
            nni, dist = connect(n, base_polyline)
            if dist < snap_radius:  # TODO criteria
                npl.append(polylines[i])
                for k in properties.keys():
                    try:
                        dummy = properties[k][i]
                    except IndexError:
                        dummy = None
                    if dummy is None:
                        print('Found in reconstruct_laterals: First order laterals that miss some properties found')
                    else:
                        nprop.setdefault(k, []).append(properties[k][i])
                        nprop["parent-poly"][-1] = 0
                        nprop["parent-node"][-1] = nni
            else:
                s += 1

    print("removed", s, "roots")
    return npl, nprop


def analyze_zones(polylines, properties):
    """
    assumes polylines[0] is tap root and determines la, lb, ln
    """
    pni = properties["parent-node"]
    ii = pni[1:]  # lateral parent node indices
    ii.sort()
#     print(ii)
    n0 = 0
    nl0 = ii[0]
    lb = polyline_length(n0, nl0, polylines[0])  # length of basal zone
    ne = ii[-1]
    nle = len(polylines[0]) - 1
    la = polyline_length(ne, nle, polylines[0])  # length of apical zone
    ln_ = []
    for i in range(0, len(ii) - 1):  # laterals
        i0 = ii[i]
        i1 = ii[i + 1]
        print(i0, i1)
        ln_.append(polyline_length(i0, i1, polylines[0]))
    return lb, ln_, la


def analyze_theta(polylines, properties):
    """
    assumes polylines[0] is tap root and determines theta
    """
#     print("analyze_theta")
#     angle = np.arccos(np.clip(np.dot(np.array([0, -1]), np.array([0, -1])), -1.0, 1.0))
#     print("down", angle)
#     angle = np.arccos(np.clip(np.dot(np.array([0, -1]), np.array([1, 0])), -1.0, 1.0))
#     print("90", angle)
#     angle = np.arccos(np.clip(np.dot(np.array([0, -1]), np.array([-1, 0])), -1.0, 1.0))
#     print("90", angle)
    theta = []
    pni = properties["parent-node"]
    for i in range(1, len(polylines)):
        j = int(pni[i])
        if j > 0:
            p = polylines[0][j]
            p0 = polylines[0][j - 1]
            if len(polylines[i]) > 2:
                p2 = polylines[i][2]
            else:
                p2 = polylines[i][0]
            v1 = p - p0
            v2 = p2 - p
            v1 = v1 / np.linalg.norm(v1)
            # v1 = [0, 0, -1] # also makes sense
            v2 = v2 / np.linalg.norm(v2)
            angle = np.arccos(np.clip(np.dot(v1, v2), -1.0, 1.0))
            theta.append(angle)
    return theta


def create_diameters(polylines, properties, functions):
    """ 
    adds diameter property as mean value of the functions
    """
    properties["diameter"] = []
    for i, p in enumerate(polylines):
        d = 0.
        for j in range(0, len(p)):
#             print(len(p))
#             print(len(functions["diameter"]))
            d += functions["diameter"][i][j] / len(p) / 10.  # mm -> cm
        properties["diameter"].append(d)


def create_age_delay(polylines, properties, time, r0, lmax0, lad0):
    """
    adds estimated lateral root age to the properties, based on given tap root paraeters 
    """
    pni = properties["parent-node"]
    properties.setdefault("age", []).append(time)  # tap root et =  0
    for i in range(1, len(polylines)):
        il = pni[i]
        bl = polyline_length(0, il, polylines[0])  # latearal base length
        et = negexp_age(bl, r0, lmax0) + lad0  # time the lateral emerged
        age = max(time - et, 0.)
        age = min(age, time)
        properties["age"].append(age)  # tap root et =  0

