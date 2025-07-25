""" coupling with DuMux as solver for the soil part, dumux-rosi must be installed & compiled """
import sys; sys.path.append("../.."); sys.path.append("../../src/")
sys.path.append("../../../dumux-rosi/build-cmake/cpp/python_binding/")  # dumux python binding
sys.path.append("../../../dumux-rosi/python/modules/")  # python wrappers

import plantbox as pb
import visualisation.vtk_plot as vp
from functional.PlantHydraulicParameters import PlantHydraulicParameters
from functional.PlantHydraulicModel import HydraulicModel_Doussan
from functional.PlantHydraulicModel import HydraulicModel_Meunier
from rosi_richards import RichardsSP  # C++ part (Dumux binding)
from richards import RichardsWrapper  # Python part
import numpy as np
import matplotlib.pyplot as plt
import timeit


def sinusoidal(t):
    return np.sin(2. * np.pi * np.array(t) - 0.5 * np.pi) + 1.


""" Parameters """  # |\label{l71c:param}|
min_b = [-35., -10., -50.]  # [cm]
max_b = [35., 10., 0.]  # [cm]
cell_number = [17, 5, 50]  # 1D

path = "../../modelparameter/structural/rootsystem/"
name = "Zeamays_synMRI_modified"  #"Anagallis_femina_Leitner_2010"  # Zea_mays_1_Leitner_2010, Zeamays_synMRI.xml  <<<<-------
trans = 500  # cm3 /day (sinusoidal) = mL/day
wilting_point = -15000  # cm
rs_age = 21  # root system initial age [day]

loam = [0.08, 0.43, 0.04, 1.6, 50]
initial = -300  # cm

sim_time = 7.5  # [day]
dt = 3600. / (24 * 3600)  # [days]  # |\label{l71c:param_end}|

""" Initialize macroscopic soil model """
s = RichardsWrapper(RichardsSP())  # |\label{l71c:soil}|
s.initialize()
s.createGrid(min_b, max_b, cell_number, periodic = True)  # [cm]
s.setHomogeneousIC(initial, True)  # [cm] total potential
s.setTopBC("noFlux")
s.setBotBC("noFlux")
s.setVGParameters([loam])
s.setParameter("Soil.SourceSlope", "1000")  # |\label{l71c:regularisation}|
s.initializeProblem()
s.setCriticalPressure(wilting_point)  # |\label{l71c:soil_end}|

""" Initialize xylem model """
plant = pb.MappedPlant()  # |\label{l71c:soil_plant}|
plant.readParameters(path + name + ".xml")
sdf = pb.SDF_PlantBox(np.inf, np.inf, max_b[2] - min_b[2] - 0.5)  # |\label{l71c:domain}|
plant.setGeometry(sdf)  # |\label{l71c:soil_plant_end}|

""" root hydraulic properties """
params = PlantHydraulicParameters()  # |\label{l71c:hydraulic}|
params.read_parameters("../../modelparameter/functional/plant_hydraulics/couvreur2012")
# params.plot_conductivities(True)
hm = HydraulicModel_Doussan(plant, params)
hm.wilting_point = wilting_point  # |\label{l71c:hydraulic_end}|

""" Coupling (map indices) """
picker = lambda x, y, z: s.pick([x, y, z])  # |\label{l71c:coupling}|
plant.setSoilGrid(picker)
plant.initialize(True)
plant.simulate(rs_age, True)
hm.test()  # |\label{l71c:test}|

""" Numerical solution """
start_time = timeit.default_timer()
t = 0.
x_, y_ = [], []
N = round(sim_time / dt)

for i in range(0, N):  # |\label{l71c:loop}|

    plant.simulate(dt)  # |\label{l71c:plant}|
    # hm.test_unmapped()
    hs = s.getSolutionHead()  # |\label{l71c:hs}|
    hx = hm.solve(rs_age + t, -trans * sinusoidal(t), hs, cells = True)  # |\label{l71c:hx}|

    fluxes = hm.soil_fluxes(rs_age + t, hx, hs)  # |\label{l71c:soil_model}|
    s.setSource(fluxes)
    s.solve(dt)  # |\label{l71c:soil_model_end}|

    x_.append(t)
    y_.append(float(hm.get_transpiration(rs_age + t, hx, hs, cells = True)))  # |\label{l71c:results}|

    n = round(float(i) / float(N) * 100.)  # |\label{l71c:progress}|
    print("[" + ''.join(["*"]) * n + ''.join([" "]) * (100 - n) + "], [{:g}, {:g}] cm soil [{:g}, {:g}] cm root at {:g} days {:g}"
            .format(np.min(hs), np.max(hs), np.min(hx), np.max(hx), s.simTime, hx[0]))

    if i % 10 == 0:  # |\label{l71c:write}|
        vp.write_soil("results/example72_{:06d}".format(i // 10), s, min_b, max_b, cell_number)
        vp.write_plant("results/example72_{:06d}".format(i // 10), hm.ms.plant())  # |\label{l71c:write_end}|

    t += dt  # [day]

print ("Coupled benchmark solved in ", timeit.default_timer() - start_time, " s")  # |\label{l71c:timing}|
# print("Len:", len(hx), len(hm.ms.mappedSegments().segments))

""" VTK visualisation """  # |\label{l71c:plots}|
vp.plot_roots_and_soil(hm.ms.mappedSegments(), "matric potential", hx, s, True, np.array(min_b), np.array(max_b), cell_number, name)

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

# kz = 4.32e-2  # [cm^3/day] # TODO have we described the options in the book?
# kr = 1.728e-4  # [1/day]
# params.set_kr_const(kr)
# params.set_kx_const(kz)
