import sys;  sys.path.append("../..")

import numpy as np
import scipy
from scipy import integrate
import matplotlib.pyplot as plt

#
# Model parameter
#
Q = 4  # µg/d/tip
Dl = 2.43e-6 * 3600 * 24  # cm2/d
theta = 0.3
R = 16.7  # 16.7  # -
k = 2.60e-6 * 3600 * 24  # d-1
L = 4  # cm (for line source only)

r = 2  # initial growth speed (linear growth)
simtime = 10  # days
depth = 26  # cm


def to32(x):
    return np.sqrt(x * x * x)


def fun11a(x, t):
    c = -R / (4 * Dl * t)
    d = 8 * theta * to32(np.pi * Dl * t)
    xtip = (-3) - (r * (simtime - t))
    z = x - xtip
    return  (Q * np.sqrt(R)) / d * np.exp(c * (z * z) - k / R * t)


def fun11b(x, t, l):
    c = -R / (4 * Dl * t)
    d = 8 * theta * to32(np.pi * Dl * t);
    xtip = (-3) - (r * (simtime - t))
    if xtip + l < 0:
        z = x - (xtip + l);
        return (Q * np.sqrt(R)) / d * np.exp(c * z * z - k / R * t)
    else:
        return 0

# print(fun11a(-3, 1e-99))
# print(fun11a(-3, 1))
# print(fun11a(-23, 10))


N = 100

# t_ = np.arange(1, 11)
# z_ = np.linspace(0, -depth, N)
# cz = np.zeros((N, 11))
# for j in t_:
#     for i, z in enumerate(z_):
#         f = lambda t: fun11a(z, t)
#         cz[i, j] = f(j)  # res[-1]
#     plt.plot(z_, cz[:, j], "g:")
# plt.ylabel("µg /d")
# plt.xlabel("z (cm)")
# plt.title("Integrand for t = 1,2,...,11 days (rootage = 10 days)")
# plt.show()
#
# t_ = np.linspace(1, 10, N)
# z_ = np.arange(1, 30, 2);
# ct = np.zeros((N, 30))
# for j in z_:# 2) Simulates root growth
#     for i, t in enumerate(t_):
#         f = lambda z: fun11a(z, t)
#         ct[i, j] = f(-j)  # res[-1]
#     plt.plot(t_, ct[:, j])
# plt.ylabel("µg /d")
# plt.xlabel("time (d)")
# plt.legend(["-1", "-3", "-5", "-7", "-9"])
# plt.title("Integrand at z = -1,-3,...,-29 cm (rootage = 10 days)")
# plt.show()

# # Moving Point Source
# z_ = np.linspace(-30, 0, N);
# c = np.zeros((len(z_),))
# for i, z in enumerate(z_):
#     f = lambda t: fun11a(z, t)
#     res = integrate.quad(f, 0.2, 10)  # <- t0 determines peak hight
#     c[i] = res[0]
#
# plt.plot(np.abs(z_), c, 'r')
# plt.plot([23], [0], 'k*') # tip
# plt.plot([3], [0], 'r*') # base
# plt.show()

# Moving Line Source
z_ = np.linspace(-30, 0, N);
c = np.zeros((len(z_),))
for i, z in enumerate(z_):
    f = lambda l, t : fun11b(z, t, l)  # care: dblquad f(y,x)
    res = integrate.dblquad(f, 0.2, 10, lambda x: 0, lambda x: L)
    c[i] = res[0]

plt.plot(np.abs(z_), c, 'r')
plt.plot([23], [0], 'k*')  # tip
plt.plot([3], [0], 'g*')  # base
plt.show()
