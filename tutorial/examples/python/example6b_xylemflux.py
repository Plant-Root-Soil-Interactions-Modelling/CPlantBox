""" water movement within the root (static soil) """
import sys; sys.path.append("../../..")
from xylem_flux import XylemFluxPython  # Python hybrid solver
import plantbox as pb
import vtk_plot as vp

import numpy as np
import matplotlib.pyplot as plt

""" Parameters """
kz = 4.32e-2  # axial conductivity [cm^3/day]
kr = 1.728e-4  # radial conductivity [1/day]
p_s = -200  # static soil pressure [cm]
p0 = -500  # dircichlet bc at top
simtime = 14  # [day] for task b

""" root system """
rs = pb.MappedRootSystem()
path = "../../../modelparameter/rootsystem/"
name = "Anagallis_femina_Leitner_2010"  # Zea_mays_1_Leitner_2010
rs.readParameters(path + name + ".xml")
rs.initialize()
rs.simulate(simtime, False)

""" root problem """
r = XylemFluxPython(rs)
r.setKr([kr])
r.setKx([kz])
nodes = r.get_nodes()
soil_index = lambda x, y, z : 0
r.rs.setSoilGrid(soil_index)

""" Numerical solution """
rx = r.solve_dirichlet(0., p0, p_s, [p_s], True)
fluxes = r.segFluxes(simtime, rx, -200 * np.ones(rx.shape), False)  # cm3/day
print("Transpiration", r.collar_flux(simtime, rx, [p_s]), "cm3/day")

""" plot results """
plt.plot(rx, nodes[:, 2] , "r*")
plt.xlabel("Xylem pressure (cm)")
plt.ylabel("Depth (m)")
plt.title("Xylem matric potential (cm)")
plt.show()

""" Additional vtk plot """
ana = pb.SegmentAnalyser(r.rs)
ana.addData("rx", rx)
ana.addData("fluxes", np.maximum(fluxes, -1.e-3))  # cut off for vizualisation
vp.plot_roots(ana, "rx", "Xylem matric potential (cm)")  # "fluxes"
