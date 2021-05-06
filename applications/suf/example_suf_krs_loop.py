""" 
soil uptake fraction of a root system (soil is in hydrostatic equilibrium) 
"""
import sys; sys.path.append("../../python/modules/"); sys.path.append("../../../CPlantBox/"); 
sys.path.append("../../../CPlantBox/src/python_modules")

from xylem_flux import XylemFluxPython  # Python hybrid solver
import plantbox as pb
import vtk_plot as vp

import numpy as np
import matplotlib.pyplot as plt


def get_age(l, r, k):
    l = np.minimum(l, k - 0.001)  # to avoid nans
    return -k / r * np.log(1 - l / k)


""" Parameters """
# artificial shoot
kr0 = np.array([[-154, 0.], [0, 1.e-12], [1e20, 1.e-12]])
kz0 = np.array([[-154, 0.], [0, 1.], [1e20, 1.]])

# Doussan et al. TODO distance to age
kr1 = np.array([[-1e4, 1.e-12], [0., 1.e-12], [0., 2.55e-6], [12.5, 2.55e-6], [20.9, 8.9e-7], [44.6, 8.9e-7], [62.7, 2.1e-7], [100, 2.1e-7]])  # negative values must be samll for kr
kr2 = np.array([[-1e4, 1.e-12], [0., 1.e-12], [0., 2.e-4], [10, 2.e-4], [15, 3.e-5], [20, 3.e-5]])
kr3 = np.array([[-1e4, 1.e-12], [0., 1.e-12], [0., 2.e-4], [10, 2.e-4], [15, 3.e-5], [20, 3.e-5]])
kz1 = np.array([[0, 2.3148e-4], [18.3, 2.3148e-4], [27.8, 4.0509e-3], [36.4, 4.0509e-3], [51.1, 5.752278e-2], [100, 5.752278e-2]])
kz2 = np.array([[0, 1.e-6], [9, 2.e-4], [13, 6.e-4], [20, 6.e-4]])
kz3 = np.array([[0, 1.e-6], [9, 2.e-4], [13, 6.e-4], [20, 6.e-4]])

# Const
# kr1 = np.array([[-154, 0.], [0, 1.728e-4], [1e20, 1.728e-4]])
# kr2 = np.array([[-154, 0.], [0, 1.728e-4], [1e20, 1.728e-4]])
# kr3 = np.array([[-154, 0.], [0, 1.728e-4], [1e20, 1.728e-4]])
# kz1 = np.array([[-154, 0.], [0, 4.32e-2], [1e20, 4.32e-2]])
# kz2 = np.array([[-154, 0.], [0, 4.32e-2], [1e20, 4.32e-2]])
# kz3 = np.array([[-154, 0.], [0, 4.32e-2], [1e20, 4.32e-2]])

# Old values ???
# kr1 = np.array([[-1e4, 1.e-9], [0., 1.e-9], [0, 1.14e-03], [2, 1.09e-03], [4, 1.03e-03], [6, 9.83e-04], [8, 9.35e-04], [10, 8.90e-04], [12, 8.47e-04], [14, 8.06e-04], [16, 7.67e-04], [18, 7.30e-04], [20, 6.95e-04], [22, 6.62e-04], [24, 6.30e-04], [26, 5.99e-04], [28, 5.70e-04], [30, 5.43e-04], [32, 5.17e-04]])
# kr2 = np.array([[-1e4, 1.e-9], [0., 1.e-9], [0, 4.11e-03], [1, 3.89e-03], [2, 3.67e-03], [3, 3.47e-03], [4, 3.28e-03], [5, 3.10e-03], [6, 2.93e-03], [7, 2.77e-03], [8, 2.62e-03], [9, 2.48e-03], [10, 2.34e-03], [11, 2.21e-03], [12, 2.09e-03], [13, 1.98e-03], [14, 1.87e-03], [15, 1.77e-03], [16, 1.67e-03], [17, 1.58e-03]])
# kr3 = np.array([[-1e4, 1.e-9], [0., 1.e-9], [0, 4.11e-03], [1, 3.89e-03], [2, 3.67e-03], [3, 3.47e-03], [4, 3.28e-03], [5, 3.10e-03], [6, 2.93e-03], [7, 2.77e-03], [8, 2.62e-03], [9, 2.48e-03], [10, 2.34e-03], [11, 2.21e-03], [12, 2.09e-03], [13, 1.98e-03], [14, 1.87e-03], [15, 1.77e-03], [16, 1.67e-03], [17, 1.58e-03]])
# kz1 = np.array([[0, 6.74e-02], [2, 7.48e-02], [4, 8.30e-02], [6, 9.21e-02], [8, 1.02e-01], [10, 1.13e-01], [12, 1.26e-01], [14, 1.40e-01], [16, 1.55e-01], [18, 1.72e-01], [20, 1.91e-01], [22, 2.12e-01], [24, 2.35e-01], [26, 2.61e-01], [28, 2.90e-01], [30, 3.21e-01], [32, 3.57e-01]])
# kz2 = np.array([[0, 4.07e-04], [1, 5.00e-04], [2, 6.15e-04], [3, 7.56e-04], [4, 9.30e-04], [5, 1.14e-03], [6, 1.41e-03], [7, 1.73e-03], [8, 2.12e-03], [9, 2.61e-03], [10, 3.21e-03], [11, 3.95e-03], [12, 4.86e-03], [13, 5.97e-03], [14, 7.34e-03], [15, 9.03e-03], [16, 1.11e-02], [17, 1.36e-02]])
# kz3 = np.array([[0, 4.07e-04], [1, 5.00e-04], [2, 6.15e-04], [3, 7.56e-04], [4, 9.30e-04], [5, 1.14e-03], [6, 1.41e-03], [7, 1.73e-03], [8, 2.12e-03], [9, 2.61e-03], [10, 3.21e-03], [11, 3.95e-03], [12, 4.86e-03], [13, 5.97e-03], [14, 7.34e-03], [15, 9.03e-03], [16, 1.11e-02], [17, 1.36e-02]])

simtime = 60  # [day] for task b

""" root system """
rs = pb.MappedRootSystem()
rs.setGeometry(pb.SDF_PlantBox(1.e6, 1.e6, 1.e6))  # not allowed to grow upwards out of soil
p_s = np.linspace(-500, -200, 3001)  #  -200.*np.ones((2001, 1))   # 3 meter down, from -200 to -500, resolution in mm
soil_index = lambda x, y, z: int(-10 * z)  # maps to p_s (hydrostatic equilibirum)
rs.setSoilGrid(soil_index)

path = "../../../CPlantBox//modelparameter/rootsystem/"
name = "Zea_mays_1_Leitner_2010"  # "Glycine_max"  # "Anagallis_femina_Leitner_2010"  # Zea_mays_1_Leitner_2010
rs.readParameters(path + name + ".xml")
rs.getRootSystemParameter().seedPos.z = -3
rs.setSeed(2)  # random

# quick fix
r_, k_ = [], []
for p in rs.getRootRandomParameter():
    p.dx = 0.1
    r_.append(p.r)
    k_.append(p.lmax)
kr1[2:, 0] = get_age(kr1[2:, 0], r_[1], k_[1]) 
kz1[:, 0] = get_age(kz1[:, 0], r_[1], k_[1]) 
kr2[2:, 0] = get_age(kr2[2:, 0], r_[2], k_[2]) 
kz2[:, 0] = get_age(kz2[:, 0], r_[2], k_[2])
kr3[2:, 0] = get_age(kr3[2:, 0], r_[3], k_[3]) 
kz3[:, 0] = get_age(kz3[:, 0], r_[3], k_[3]) 

rs.initialize()
rs.simulate(simtime, False)

""" set up xylem parameters """
r = XylemFluxPython(rs)
# defaults...
kr4 = kr1  # basal
kr5 = kr1  # shoot borne
kz4 = kz1
kz5 = kz1
r.setKrTables([kr0[:, 1], kr1[:, 1], kr2[:, 1], kr3[:, 1], kr4[:, 1], kr5[:, 1]], [kr0[:, 0], kr1[:, 0], kr2[:, 0], kr3[:, 0], kr4[:, 0], kr5[:, 0]])
r.setKxTables([kz0[:, 1], kz1[:, 1], kz2[:, 1], kz3[:, 1], kz4[:, 1], kz5[:, 1]], [kz0[:, 0], kz1[:, 0], kz2[:, 0], kz3[:, 0], kz4[:, 0], kz5[:, 0]])

""" for debugging """
r.test()
r.plot_conductivities()
shoot_segs = rs.getShootSegments()
print("Shoot segments", [str(s) for s in shoot_segs])
print("Shoot type", rs.subTypes[0])

jc_, krs_, l_, eswp_, suf_ = [], [], [], [], []
""" numerical solution of transpiration -1 cm3/day"""
simtime += 1
for t in range(10, simtime):
          
    suf = r.get_suf(t)
    krs, _ = r.get_krs(t)    
    krs_.append(krs)
    
    if t > 10 and t % 10 == 0:
        print(t)
        suf_.append(suf)

time = np.linspace(10, simtime, simtime - 10)
plt.plot(time, krs_, label="krs")
plt.xlabel("simulation time")
plt.ylabel("root system conductance $krs$ $[cm^2 d^{-1}]$")
plt.legend()
plt.show()

""" SUF plot """
fig, (ax) = plt.subplots(1, len(suf_))
ax[0].set_ylabel("depth [cm]")
for i, suf in enumerate(suf_):
    cols = ["r", "g", "b", "m"]
    x_ = np.linspace(0, -100, 100)    
    for j in range(0, 4):
        ana = pb.SegmentAnalyser(r.rs) 
        ana.addData("SUF", suf)    
        ana.filter("subType", j + 1)
        y_ = ana.distribution("SUF", 0, -100, 100, True)
        ax[i].plot(y_, x_, cols[j])
    ax[i].set_title(str(20 + 10 * i))
    ax[i].set_xlabel("SUF [-]")
    ax[i].legend([str(j) for j in range(1, 5)])
plt.show()

