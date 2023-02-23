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

""" Parameters """
# artificial shoot
kr0 = np.array([[0, 1.e-20], [1e20, 1.e-20]]) 		# kr of tap root at age 0
kz0 = np.array([[0, 0.356832], [1e20, 0.356832]])		# kx of tap root at age 1e20

# tap root
kr1 = np.array([[0, 1.14048e-3], [2, 1.08864e-3], [4, 1.0368e-3], [6, 9.8486e-4], [8, 9.3312e-4], [10, 8.8992e-4], [12, 8.47584e-4], [14, 8.06112e-4], [16, 7.67232e-4], [18, 7.3008e-4], [20, 6.9552e-4], [22, 6.61824e-4], [24, 6.29856e-4], [26, 5.99616e-4], [28, 5.7024e-4], [30, 5.42592e-4], [32, 5.16672e-4], [1e20, 5.16672e-4]])
kz1 = np.array([[0, 0.067392], [2, 0.074736], [4, 0.082944], [6, 0.092448], [8, 0.101952], [10, 0.113184], [12, 0.126144], [14, 0.139968], [16, 0.154656], [18, 0.171936], [20, 0.190944], [22, 0.21168], [24, 0.235008], [26, 0.260928], [28, 0.28944], [30, 0.321408], [32, 0.356832], [1e20, 0.356832]])

# first order lateral
kr2 = np.array([[0, 4.11264e-3], [1, 3.888e-3], [2, 3.672e-3], [3, 3.47328e-3], [4, 3.2832e-3], [5, 3.10176e-3], [6, 2.92896e-3], [7, 2.77344e-3], [8, 2.61792e-3], [9, 2.47968e-3], [10, 2.34144e-3], [11, 2.21184e-3], [12, 2.09952e-3], [13, 1.97856e-3], [14, 1.86624e-3], [15, 1.76256e-3], [16, 1.66752e-3], [17, 1.58112e-3], [1e20, 1.58112e-3]])
kz2 = np.array([[0, 4.06944e-4], [1, 5.00256e-4], [2, 6.15168e-4], [3, 7.56e-4], [4, 9.3312e-4], [5, 1.14048e-3], [6, 1.40832e-3], [7, 1.728e-3], [8, 2.12544e-3], [9, 2.60928e-3], [10, 3.21408e-3], [11, 3.94848e-3], [12, 4.85568e-3], [13, 5.97024e-3], [14, 7.344e-3], [15, 8.9856e-3], [16, 0.0110592], [17, 0.0136512], [1e20, 0.0136512]])

# second order lateral
kr3 = np.array([[0, 4.11264e-3], [1, 3.888e-3], [2, 3.672e-3], [3, 3.47328e-3], [4, 3.2832e-3], [5, 3.10176e-3], [6, 2.92896e-3], [7, 2.77344e-3], [8, 2.61792e-3], [9, 2.47968e-3], [10, 2.34144e-3], [11, 2.21184e-3], [12, 2.09952e-3], [13, 1.97856e-3], [14, 1.86624e-3], [15, 1.76256e-3], [16, 1.66752e-3], [17, 1.58112e-3], [1e20, 1.58112e-3]])
kz3 = np.array([[0, 4.06944e-4], [1, 5.00256e-4], [2, 6.15168e-4], [3, 7.56e-4], [4, 9.3312e-4], [5, 1.14048e-3], [6, 1.40832e-3], [7, 1.728e-3], [8, 2.12544e-3], [9, 2.60928e-3], [10, 3.21408e-3], [11, 3.94848e-3], [12, 4.85568e-3], [13, 5.97024e-3], [14, 7.344e-3], [15, 8.9856e-3], [16, 0.0110592], [17, 0.0136512], [1e20, 0.0136512]])

simtime = 154  # [day] for task b

""" root system """
rs = pb.MappedRootSystem()
rs.setGeometry(pb.SDF_PlantBox(1.e6, 1.e6, 1.e6))  # not allowed to grow upwards out of soil
p_s = np.linspace(-500, -200, 3001)  #  -200.*np.ones((2001, 1))   # 3 meter down, from -200 to -500, resolution in mm
soil_index = lambda x, y, z: int(-10 * z)  # maps to p_s (hydrostatic equilibirum)
rs.setSoilGrid(soil_index)

path = "../../../CPlantBox//modelparameter/rootsystem/"
name = "Glycine_max_Moraes2020_opt2"  # "Glycine_max"  # "Anagallis_femina_Leitner_2010"  # Zea_mays_1_Leitner_2010
rs.readParameters(path + name + ".xml")
rs.getRootSystemParameter().seedPos.z = -3
rs.setSeed(2)
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

krs_, suf_ = [], []
""" numerical solution of transpiration -1 cm3/day"""
simtime += 1
for t in range(10, simtime):
          
    suf = r.get_suf(t)
    krs, _ = r.get_krs(t)    
    krs_.append(krs)
    
    if t > 10 and t % 30 == 0:
        print(t)
        suf_.append(suf)

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
    ax[i].set_title(str(30 + 30 * i))
    ax[i].set_xlabel("SUF [-]")
    ax[i].legend([str(j) for j in range(1, 5)])
plt.savefig("results/" + name + "/" + name + "_SUF.pdf", dpi = 300, bbox_inches='tight')
plt.show()

