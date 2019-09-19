#
# Length Benchmark
#
# Compares the analytically calculated lengths with the rootbox approximation
# parameters are set within the code L65- (with no standard deviation)
#
# Benchmark 1: single root, no laterals
# Benchmark 2: single root with laterals
# Benchmark 3: basal roots, no laterals
# Benchmark 4: basal roots with laterals
#
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

import ../rootbox as rb


def maxRootLength(la, lb, ln, nob):  # maximal length the root will reach
    return la + lb + (nob - 1) * ln


def rootLength(t, r, k):  # root length at a certain age
    return k * (1 - np.exp(-r * t / k))


def rootAge(l, r, k):  # root age at a certain length
    return -np.log(1 - l / k) * k / r


def rootLateralLength(t, et, r, k):  # length of first order laterals (without second order laterals)
    i = 0
    l = 0
    while et[i] < t:
        age = t - et[i]
        l += rootLength(age, r, k)
        i += 1
    return l


def v2a(vd):  # rb.std_vector_double_ to numpy array
    l = np.zeros(len(vd))
    for i in range(0, len(vd)):
        l[i] = vd[i]
    return l


def a2v(a):  #  numpy array to rb.std_vector_double
    l = rb.std_vector_double_()
    for d in a:
        l.append(d)
    return l


def a2i(a):  #  numpy array to rb.std_vector_int
    l = rb.std_vector_int_()
    for i in a:
        l.append(i)
    return l


def vv2a(vd):  # rb.std_vector_Vector3_ to numpy array
    N = len(vd)
    l = np.zeros((N, 3))
    for i in range(0, N):
        l[i, :] = [vd[i].x, vd[i].y, vd[i].z]
    return l


#
# Root type parameter
#
rs = rb.RootSystem()
p0 = rb.RootTypeParameter(rs)
p1 = rb.RootTypeParameter(rs)

# Taproot
p0.name = "taproot"
p0.type = 1
p0.lb = 1
p0.la = 10
p0.nob = 20
p0.ln = 89. / 19.
p0.r = 1
p0.dx = 0.5
p0.k = maxRootLength(p0.la, p0.lb, p0.ln, p0.nob)
# print(p0)

# 1st order lateral
p1.name = "lateral"
p1.type = 2
p1.la = 25
p1.ln = 0
p1.r = 2
p1.k = 25
p1.dx = 0.1
# print(p1)

#
# Root system parameter (neglecting shoot borne)
#
maxB = 100
firstB = 10.
delayB = 3.
rsp = rb.RootSystemParameter()
rsp.set(-3., firstB, delayB, maxB, 0, 1.e9, 1.e9, 1.e9, 0., 0.)

times = np.array([0., 7., 15., 30., 60.])
dt = np.diff(times)
times = times[1:]

#
# Benchmark 1: single root, no laterals
#
print("* ")
print("* Benchmark 1: single root, no laterals")
print("*")

# Analytical
l = rootLength(times, p0.r, p0.k)

# Numerical
rs.setOrganTypeParameter(p0)
rs.setOrganTypeParameter(p1)
rs.initialize()
c = 0
nl = np.zeros(len(times))
non = np.zeros(len(times))
for t in dt:
    rs.simulate(t, False)
    d = v2a(rs.getScalar(rb.ScalarType.length))
    nl[c] = d[0]
    non[c] = rs.getNumberOfNodes()
    c += 1

# Numerical, same, but with tiny time stepping
rs.initialize()  # resets everything
c = 0
nl2 = np.zeros(len(times))
non2 = np.zeros(len(times))
for t in dt:
    dt_ = t / 1000.
    for i in range(0, 1000):
        rs.simulate(dt_, False)
    d = v2a(rs.getScalar(rb.ScalarType.length))
    nl2[c] = d[0]
    non2[c] = rs.getNumberOfNodes()
    c += 1

print("times \t\t\t", times)
print("analytical lenghts \t", l)
print("numerical lenghts \t", nl)
print("numerical lenghts 2 \t", nl2, "\n")

print("mean axial resoltuion = \t", nl / non)  # should be less but in the same order as the defined axial resolution
print("mean axial resoltuion 2 = \t", nl2 / non2, "\n")  # should be less but in the same order as the defined axial resolution

#
# Benchmark 2: single root
#
print("* ")
print("* Benchmark 2: single root")
print("*")

# Analytical
i = 0
et = np.zeros(int(p0.nob))
while i < p0.nob:
    et[i] = rootAge(p0.la + p0.lb + p0.ln * i, p0.r, p0.k + 1e-12)
    i += 1
# print("lateral emergence times", et)

l = rootLength(times, p0.r, p0.k)  # zero order lengths
l1 = np.zeros(times.size)
j = 0
for t in times :
    l1[j] = rootLateralLength(t, et, p1.r, p1.k)
    j = j + 1

# Numerical
p0.successor = a2i([2])  # add successors
p0.successorP = a2v([1])
rs.initialize()
c = 0
nl = np.zeros(len(times))
nl0 = np.zeros(len(times))
non = np.zeros(len(times))
for t in dt:
    rs.simulate(t, False)
    d = v2a(rs.getScalar(rb.ScalarType.length))
    nl[c] = sum(d)
    nl0[c] = d[0]  # first entry is the tap root
    non[c] = rs.getNumberOfNodes()
    c += 1

# Numerical, same, but with tiny time stepping
rs.initialize()  # resets everything
c = 0
nl2 = np.zeros(len(times))
nl02 = np.zeros(len(times))
non2 = np.zeros(len(times))
for t in dt:
    dt_ = t / 1000.
    for i in range(0, 1000):
        rs.simulate(dt_, False)
    d = v2a(rs.getScalar(rb.ScalarType.length))
    nl2[c] = sum(d)
    nl02[c] = d[0]  # first entry is the tap root
    non2[c] = rs.getNumberOfNodes()
    c += 1

print("times \t\t\t\t", times)
print("analytical zero order length \t", l)
print("numerical zero order length \t", nl0)
print("analytical first order length \t", l1)
print("numerical first order length \t", nl - nl0)
print("analytical total length \t", l + l1)
print("numerical total length \t\t", nl)
print("numerical total length 2 \t", nl2, "\n")

print("mean axial resoltuion = \t", nl / non)  # should be less but in the same order as the defined axial resolution
print("mean axial resoltuion 2 = \t", nl2 / non2, "\n")  # should be less but in the same order as the defined axial resolution

#
# Benchmark 3: basal roots, no laterals
#
print("* ")
print("* Benchmark 3: basal roots, no laterals")
print("* ")

# Analytical
etB = np.array(range(maxB)) * delayB + np.ones(maxB) * firstB  # basal root emergence times
bl = np.zeros(times.size)
j = 0  # time counter
for t in times:
    i = 0  # basal root counter
    while t - etB[i] > 0:
        bl[j] += rootLength(t - etB[i], p0.r, p0.k)
        i += 1
    j += 1

# Numerical
p0.successor = a2i([])  # remove successors
p0.successorP = a2v([])
rs.setRootSystemParameter(rsp)
rs.initialize()
c = 0
nl = np.zeros(len(times))
non = np.zeros(len(times))
nl_tap = np.zeros(len(times))
nl_basal = np.zeros(len(times))
for t in dt:
    rs.simulate(t, False)
    d = v2a(rs.getScalar(rb.ScalarType.length))
    nl[c] = sum(d)
    non[c] = rs.getNumberOfNodes()

    seg = rs.getSegments()
    nodes = rs.getNodes()

    ana = rb.SegmentAnalyser(rs)
    ana.filter("subType", 1.)  # 1 is the type number of the tap root
    nl_tap[c] = ana.getSummed("length")

    ana = rb.SegmentAnalyser(rs)
    ana.filter("subType", 4.)  # 4 is the default type number of basal roots
    nl_basal[c] = ana.getSummed("length")

    c += 1

print("times \t\t\t\t", times)
print("analytical tap root lenght \t", l)
print("numerical tap root lenght \t", nl_tap)
print("analytical summed basal length \t", bl)
print("numerical summed basal length \t", nl_basal)
print("analytical total length \t", l + bl)
print("numerical total length \t\t", nl_tap + nl_basal, " (SegmentAnalyser) \n \t\t\t\t", nl, " (RootSystem)\n")

print("mean axial resoltuion = \t", nl / non)  # should be less but in the same order as the defined axial resolution

#
# Benchmark 4: basal roots, with laterals
#
print("\n* ")
print("* Benchmark 4")
print("* ")

# Analytical
# etB as berfore
# et a before
bl = np.zeros(times.size)
j = 0  # time counter
for t in times:
    i = 0  # basal root counter
    while t - etB[i] > 0:
        bl[j] += (rootLateralLength(t - etB[i], et, p1.r, p1.k) + rootLength(t - etB[i], p0.r, p0.k))
        i += 1
    j += 1

# Numerical
p0.successor = a2i([2])  # add successors
p0.successorP = a2v([1])
rs.initialize()

nl = np.zeros(len(times))
non = np.zeros(len(times))
nl_tap = np.zeros(len(times))
nl_taplateral = np.zeros(len(times))
nl_basal = np.zeros(len(times))
nl_basallateral = np.zeros(len(times))
for c, t in enumerate(dt):
    rs.simulate(t, False)
    d = v2a(rs.getScalar(rb.ScalarType.length))
    nl[c] = sum(d)
    non[c] = rs.getNumberOfNodes()

    print("Base root sub type ", rs.getBaseRoots()[1].param().subType)
    ana = rb.SegmentAnalyser(rs)
    pp = ana.getParameter("subType")
#     for p in pp:
#         print(p)

    print("one: ", ana.getSummed("one"))
    ana.filter("subType", 1.)  # 1 is the type number of the tap root
    nl_tap[c] = ana.getSummed("length")
    print("subType1: ", ana.getSummed("one"))

    ana = rb.SegmentAnalyser(rs)
    ana.filter("parentType", 1.)  # 1 is the type number of the tap root
    nl_taplateral[c] = ana.getSummed("length")

    ana = rb.SegmentAnalyser(rs)
    ana.filter("subType", 4.)  # 4 is the default type number of basal roots
    nl_basal[c] = ana.getSummed("length")

    ana = rb.SegmentAnalyser(rs)
    ana.filter("parentType", 4.)  # 4 is the default type number of basal roots
    nl_basallateral[c] = ana.getSummed("length")

print("times \t\t\t\t", times)
print("analytical tap root lenght \t", l + l1, " (...all including their laterals) ")
print("numerical tap root lenght \t", nl_tap + nl_taplateral)
print("analytical summed basal length\t", bl)
print("numerical summed basal length \t", nl_basal + nl_basallateral)
print("analytical total length \t", l + l1 + bl)
print("numerical total length  \t", nl, "\n")

print("mean axial resoltuion = \t", nl / non)  # should be less but in the same order as the defined axial resolution

# #
# # Matlab like plotting of analytical sotlution of a single root with laterals
# #
# t_ = np.linspace(0,100,100)
# l_ = rootLength(t_,p0.r,p0.k)
# l1_ = np.zeros(t_.size)
# i = 0
# for t in t_:
#     l1_[i] = rootLateralLength(t, et, p1.r, p1.k)
#     i += 1
# plt.plot(t_,l_)
# plt.plot(t_,l1_,'g')
# plt.plot(times,l,'ro')
# h0 = mpatches.Patch(color='blue', label='zero order')
# h1 = mpatches.Patch(color='green', label='first order')
# plt.xlabel("Age (days)")
# plt.ylabel("Length (cm)")
# plt.legend([h0, h1],['zero order','first order'])
# plt.show()

