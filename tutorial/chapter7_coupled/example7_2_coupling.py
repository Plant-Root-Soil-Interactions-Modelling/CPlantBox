""" coupling with DuMux as solver for the soil part, dumux-rosi must be installed & compiled """
import sys; sys.path.append("../.."); sys.path.append("../../src/")
sys.path.append("../../../dumux-rosi/build-cmake/cpp/python_binding/")  # dumux python binding
sys.path.append("../../../dumux-rosi/python/modules/")  # python wrappers

import plantbox as pb
import visualisation.vtk_plot as vp
from functional.PlantHydraulicParameters import PlantHydraulicParameters  # |\label{l41:imports}|
from functional.PlantHydraulicModel import HydraulicModel_Doussan
from functional.PlantHydraulicModel import HydraulicModel_Meunier  # |\label{l41:imports_end}|
from rosi_richards import RichardsSP  # C++ part (Dumux binding)
from richards import RichardsWrapper  # Python part
import numpy as np
import matplotlib.pyplot as plt
import timeit


def sinusoidal(t):
    return np.sin(2. * np.pi * np.array(t) - 0.5 * np.pi) + 1.


""" Parameters """  # |\label{l71c:param}|
# min_b = [-38., -8., -50.]  # [cm]
# max_b = [38., 8., 0.]  # [cm]
# # cell_number = [1, 1, 50]  # [1] spatial resolution (1D model)
# cell_number = [19, 4, 50]  # 1D

min_b = [-4., -4., -25.]
max_b = [4., 4., 0.]
cell_number = [8, 8, 25]  # [16, 16, 30]  # [32, 32, 60]

path = "../../modelparameter/structural/rootsystem/"
name = "Anagallis_femina_Leitner_2010"  #"Anagallis_femina_Leitner_2010"  # Zea_mays_1_Leitner_2010, Anagallis_femina_Leitner_2010 TODO <<<<-------
loam = [0.08, 0.43, 0.04, 1.6, 50]
initial = -600  # cm

trans = 6.4  # cm3 /day (sinusoidal) = mL/day
wilting_point = -10000  # cm

sim_time = 7.5  # [day]
rs_age = 14  # root system initial age [day]
dt = 360. / (24 * 3600)  # [days] Time step must be very small # |\label{l71c:param_end}|

""" Initialize macroscopic soil model """
s = RichardsWrapper(RichardsSP())  # |\label{l71c:soil}|
s.initialize()
s.createGrid(min_b, max_b, cell_number, True)  # [cm]
s.setHomogeneousIC(initial, True)  # [cm] total potential
s.setTopBC("noFlux")
s.setBotBC("noFlux")
s.setVGParameters([loam])
# s.setParameter("Newton.EnableChop", "True")
# s.setParameter("Newton.EnableAbsoluteResidualCriterion", "True")
s.setParameter("Soil.SourceSlope", "200")  # |\label{l71c:regularisation}|
s.initializeProblem()
s.setCriticalPressure(wilting_point)  # |\label{l71c:soil_end}|

""" Initialize xylem model """
plant = pb.MappedPlant()  # |\label{l71c:soil_plant}|
plant.readParameters(path + name + ".xml")
sdf = pb.SDF_PlantBox(np.inf, np.inf, max_b[2] - min_b[2] - 0.1)  # |\label{l71c:domain}|
plant.setGeometry(sdf)  # |\label{l71c:soil_plant_end}|

""" root hydraulic properties """
params = PlantHydraulicParameters()  # |\label{l71c:hydraulicparams}|
# params.read_parameters("../../modelparameter/functional/plant_hydraulics/couvreur2012")
kz = 4.32e-2  # [cm^3/day]
kr = 1.728e-4  # [1/day]
params.set_kr_const(kr)
params.set_kx_const(kz)
# params.plot_conductivities(True)  # for maize
hm = HydraulicModel_Meunier(plant, params)  # |\label{l71c:hydraulicparams_end}|
hm.wilting_point = wilting_point

""" Coupling (map indices) """
picker = lambda x, y, z: s.pick([x, y, z])  # |\label{l71c:coupling}|
plant.setSoilGrid(picker)  # maps segments
plant.initialize()  # |\label{l71c:init1}|
plant.simulate(rs_age, False)  # |\label{l71c:init2}|
hm.test()  # sanity checks # |\label{l71c:test}|

""" Numerical solution """
start_time = timeit.default_timer()  # |\label{l71c:init3}|
x_, y_ = [], []
sx = s.getSolutionHead()  # # |\label{l71c:ic}|
N = round(sim_time / dt)
t = 0.

for i in range(0, N):  # |\label{l71c:loop}|

    # plant.simulate(dt) # |\label{l71c:plant}|
    sx = s.getSolutionHead()  # |\label{l71c:soil}|
    rx = hm.solve(rs_age + t, -trans * sinusoidal(t), sx, cells = True)  # |\label{l71c:xylem}|

    fluxes = hm.soil_fluxes(rs_age + t, rx, sx)  # |\label{l71c:soil_model}|
    s.setSource(fluxes)
    s.solve(dt)  # |\label{l71c:soil_model_end}|

    x_.append(t)
    y_.append(float(hm.get_transpiration(rs_age + t, rx, sx, cells = True)))  # |\label{l71c:results}|

    n = round(float(i) / float(N) * 100.)  # |\label{l71c:progress}|
    print("[" + ''.join(["*"]) * n + ''.join([" "]) * (100 - n) + "], [{:g}, {:g}] cm soil [{:g}, {:g}] cm root at {:g} days {:g}"
            .format(np.min(sx), np.max(sx), np.min(rx), np.max(rx), s.simTime, rx[0]))

    t += dt  # [day]

print ("Coupled benchmark solved in ", timeit.default_timer() - start_time, " s")  # |\label{l71c:timing}|

""" VTK visualisation """  # |\label{l71c:plots}|
vp.plot_roots_and_soil(hm.ms, "pressure head", rx, s, True, np.array(min_b), np.array(max_b), cell_number, name)

""" Transpiration over time """
fig, ax1 = plt.subplots()
ax1.plot(x_, trans * sinusoidal(x_), 'k')  # potential transpiration
ax1.plot(x_, -np.array(y_), 'g')  # actual transpiration (neumann)
ax2 = ax1.twinx()
ax2.plot(x_, np.cumsum(-np.array(y_) * dt), 'c--')  # cumulative transpiration (neumann)
ax1.set_xlabel("Time [d]")
ax1.set_ylabel("Transpiration $[mL d^{-1}]$ per plant")
ax1.legend(['Potential', 'Actual', 'Cumulative'], loc = 'upper left')
np.savetxt(name, np.vstack((x_, -np.array(y_))), delimiter = ';')
plt.show()

