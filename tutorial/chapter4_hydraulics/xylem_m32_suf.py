import sys; sys.path.append("../modules"); sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src")

import plantbox as pb

from functional.xylem_flux import XylemFluxPython  # Python hybrid solver
from functional.PlantHydraulicModel import HydraulicModel_Doussan
from functional.PlantHydraulicModel import HydraulicModel_Meunier
from functional.PlantHydraulicParameters import PlantHydraulicParameters
from functional.Perirhizal import PerirhizalPython as Perirhizal

import visualisation.vtk_plot as vp

import matplotlib.pyplot as plt
import numpy as np

""" 
Benchmark M3.2 Root system: SUF and Krs
"""
min_b = [-4., -4., -15.]
max_b = [4., 4., 0.]
cell_number = [1, 1, 30]


def suf_krs(r, param):

    r.ms.setRectangularGrid(pb.Vector3d(min_b[0], min_b[1], min_b[2]),
        pb.Vector3d(max_b[0], max_b[1], max_b[2]),
        pb.Vector3d(cell_number[0], cell_number[1], cell_number[2]), False)

    simtime = 14  # [day] for task b

    """ SUF, Krs (a) """
    kr = 1.728e-4  # [1/day]
    kz = 4.32e-2  # [cm^3/day]
    param.set_kr_const(kr)
    param.set_kx_const(kz)
    suf_a = r.get_suf(simtime)
    krs_a, _ = r.get_krs(simtime)

    """ SUF, Krs (b) """
    kr0 = np.array([[0, 1.14e-03], [2, 1.09e-03], [4, 1.03e-03], [6, 9.83e-04], [8, 9.35e-04], [10, 8.90e-04], [12, 8.47e-04], [14, 8.06e-04], [16, 7.67e-04], [18, 7.30e-04], [20, 6.95e-04], [22, 6.62e-04], [24, 6.30e-04], [26, 5.99e-04], [28, 5.70e-04], [30, 5.43e-04], [32, 5.17e-04]])
    kr1 = np.array([[0, 4.11e-03], [1, 3.89e-03], [2, 3.67e-03], [3, 3.47e-03], [4, 3.28e-03], [5, 3.10e-03], [6, 2.93e-03], [7, 2.77e-03], [8, 2.62e-03], [9, 2.48e-03], [10, 2.34e-03], [11, 2.21e-03], [12, 2.09e-03], [13, 1.98e-03], [14, 1.87e-03], [15, 1.77e-03], [16, 1.67e-03], [17, 1.58e-03]])
    param.set_kr_age(kr0[:, 0], kr0[:, 1], subType = 1)
    param.set_kr_age(kr1[:, 0], kr1[:, 1], subType = [2, 3])
    kx0 = np.array([[0, 6.74e-02], [2, 7.48e-02], [4, 8.30e-02], [6, 9.21e-02], [8, 1.02e-01], [10, 1.13e-01], [12, 1.26e-01], [14, 1.40e-01], [16, 1.55e-01], [18, 1.72e-01], [20, 1.91e-01], [22, 2.12e-01], [24, 2.35e-01], [26, 2.61e-01], [28, 2.90e-01], [30, 3.21e-01], [32, 3.57e-01]])
    kx1 = np.array([[0, 4.07e-04], [1, 5.00e-04], [2, 6.15e-04], [3, 7.56e-04], [4, 9.30e-04], [5, 1.14e-03], [6, 1.41e-03], [7, 1.73e-03], [8, 2.12e-03], [9, 2.61e-03], [10, 3.21e-03], [11, 3.95e-03], [12, 4.86e-03], [13, 5.97e-03], [14, 7.34e-03], [15, 9.03e-03], [16, 1.11e-02], [17, 1.36e-02]])
    param.set_kx_age(kx0[:, 0], kx0[:, 1], subType = 1)
    param.set_kx_age(kx1[:, 0], kx1[:, 1], subType = [2, 3])
    suf_b = r.get_suf(simtime)
    krs_b, _ = r.get_krs(simtime)

    return krs_a, suf_a, krs_b, suf_b


" Solve for both models"
param = PlantHydraulicParameters()

r = HydraulicModel_Doussan("../../grids/RootSystem.rsml", param, cached = True)  # or HydraulicModel_Doussan, HydraulicModel_Meunier
krs_ad, suf_ad, krs_bd, suf_bd = suf_krs(r, param)
r = HydraulicModel_Meunier("../../grids/RootSystem.rsml", param, cached = True)  # or HydraulicModel_Doussan, HydraulicModel_Meunier
krs_am, suf_am, krs_bm, suf_bm = suf_krs(r, param)

print("\nDoussan")
print("Sum of SUF", np.sum(suf_ad), "from", np.min(suf_ad), "to", np.max(suf_ad), "summed positive", np.sum(suf_ad[suf_ad >= 0]))
print("Krs", krs_ad)
print("Sum of SUF", np.sum(suf_bd), "from", np.min(suf_bd), "to", np.max(suf_bd), "summed positive", np.sum(suf_bd[suf_bd >= 0]))
print("Krs", krs_bd)

print("\nMeunier")
print("Sum of SUF", np.sum(suf_am), "from", np.min(suf_am), "to", np.max(suf_am), "summed positive", np.sum(suf_am[suf_am >= 0]))
print("Krs", krs_am)
print("Sum of SUF", np.sum(suf_bm), "from", np.min(suf_bm), "to", np.max(suf_bm), "summed positive", np.sum(suf_bm[suf_bm >= 0]))
print("Krs", krs_bm)

print("cell -0.5 cm", r.ms.soil_index(0., 0., -0.5))  # first cell is at position 14
print("cell -14.5 cm", r.ms.soil_index(0, 0, -14.5))

fig, (ax1, ax2) = plt.subplots(1, 2)

z_ = np.linspace(-15, 0, cell_number[2])
z_seg = r.ms.getSegmentZ()

peri = Perirhizal(r.ms)
suf_ad_layers = np.array(peri.aggregate(suf_ad))
suf_am_layers = np.array(peri.aggregate(suf_am))
suf_bd_layers = np.array(peri.aggregate(suf_bd))
suf_bm_layers = np.array(peri.aggregate(suf_bm))

""" SUF per layer """
ax1.plot(suf_ad_layers, z_, "*", label = "(a) doussan")
ax1.plot(suf_am_layers, z_, "*", label = "(a) meunier")
ax2.plot(suf_bd_layers, z_, "*", label = "(b) doussan")
ax2.plot(suf_bm_layers, z_, "*", label = "(b) meunier")
""" SUF per segment """
# ax1.plot(suf_ad, z_seg, "*")
# ax1.plot(suf_am, z_seg, "*")
# ax2.plot(suf_bd, z_seg, "*")
# ax2.plot(suf_bm, z_seg, "*")
ax1.set_ylabel("depth [cm]")
ax1.set_xlabel("SUF [1]")
ax2.set_xlabel("SUF [1]")

ax1.legend()
ax2.legend()
plt.show()
