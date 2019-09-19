import numpy as np
import scipy.sparse.linalg as LA
from scipy import sparse
import matplotlib.pylab as plt

import xylem_flux
import timeit

from math import *

g = 9.8  # gravitational acceleration (m/s^2)
rho = 1.e3  # density of water, (kg/m^3)
ref = 0


def toPa(ph):  # cm pressure head to Pascal (kg/ (m s^2))
    return ref - ph / 100. * rho * g


def toHead(pa):  # Pascal (kg/ (m s^2)) to cm pressure head
    return -(pa - ref) * 100 / rho / g

#
# Compares analytic and numerical solution of a single root
#

# Parameters


L = 0.5  # length of single straight root (m)
a = 2.e-3  # radius (m)
kz = 5.e-13  # axial conductivity (m^5 s / kg) (mal rho ergibt die alten einheiten)
kr = 2.e-9  # radial conductivity per root type (m^2 s / kg)

p0 = toPa(-1000)  # dircichlet bc at top (ćm)
pL = toPa(-500)  # dircichlet bc at bot (ćm)
p_s = toPa(-200)  # static soil pressure (cm)

#
# analytical solution
#
# the solution of
# d^2/dz^2 p_r = - c p_r + c p_s,
# is
# p_r = p_s + d_1 exp(sqrt(c) z ) + d_2 exp(-sqrt(c) z)
# with
# c = 2 a pi kr/kz
#

c = 2 * a * pi * kr / kz

# BC
# top: p_r(0) = p0
# bot: p_r(-L) = pL
# bot: qz(L) = 0, -> d/dz p_r (L) = rho*g

AA = np.array([[1, 1], [sqrt(c) * exp(-sqrt(c) * L), -sqrt(c) * exp(sqrt(c) * L)] ])  # dirichlet top, neumann bot
bb = np.array([p0 - p_s, -rho * g])  #

# AA = np.array([[1,1], [exp(-sqrt(c)*L), exp(sqrt(c)*L)] ]) # dirichlet top & bot
# bb = np.array([p0-p_s, pL-p_s]) #

d = np.linalg.solve(AA, bb)  # compute constants d_1 and d_2 from bc

p_r = lambda z: toHead(p_s + d[0] * exp(sqrt(c) * z) + d[1] * exp(-sqrt(c) * z))

za_ = np.linspace(0, -L, 100)
pr = list(map(p_r, za_))
print(0, p_r(0), pr[0])
print(L, p_r(-L), pr[-1])

print(toHead(d))
print(c)

#
# numerical solutionprint("fin")
#

nnz = 100

# create grid
nodes = np.zeros((nnz, 3))
seg = np.zeros(((nnz - 1), 2), dtype = int)
z_ = np.zeros((nnz, 1))
c = 0
for i in range(1, nnz):
    seg[c, 0] = i - 1
    seg[c, 1] = i
    c += 1
    nodes[i, :] = [0., 0., -i * L / (nnz - 1)]
    z_[i] = -i * L / (nnz - 1)

# from constant to per segment
kr_ = [kr] * 2 * (nnz - 1)
kz_ = [kz] * 2 * (nnz - 1)
a_ = [a] * 2 * (nnz - 1)

# call back function for soil potential
soil = lambda x, y, z : p_s

# calculate fluxes within the root system
Q, b = xylem_flux.linear_system(seg, nodes, a_, kr_, kz_, rho, g, soil)  #

# Q, b = xylem_flux.bc_dirichlet(Q, b, np.array([0,nnz-1]), np.array([p0,pL]))
Q, b = xylem_flux.bc_dirichlet(Q, b, np.array([0]), np.array([p0]))  # dirichlet top
Q, b = xylem_flux.bc_neumann(Q, b, [nnz - 1], [0])  # neumann tip
# plt.spy(Q)
# plt.show()

start = timeit.default_timer()
x = LA.spsolve(Q, b, use_umfpack = True)  # direct
stop = timeit.default_timer()
# print ("linear system solved in", stop - start, " s")

# plot results
plt.plot(list(map(toHead, x)), z_, "r*")
plt.plot(pr, za_)
plt.xlabel("Xylem pressure (cm)")
plt.ylabel("Depth (m)")
plt.show()

print("done.")

