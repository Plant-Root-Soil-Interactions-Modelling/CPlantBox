import sys; sys.path.append("../../.."); sys.path.append("../../../src/python_modules")
import plantbox as pb
import vtk_plot as vp
from xylem_flux import XylemFluxPython  # Python hybrid solver
import numpy as np
import matplotlib.pyplot as plt
from Leuning import Leuning

simtime = 14  # [day] 
pl = pb.MappedPlant() #for plant objects
path = "../../../modelparameter/plant/" 
name = "manyleaves" 
pl.readParameters(path + name + ".xml")

soilSpace = pb.SDF_PlantContainer(500,500,36, True) #to avoid root growing aboveground
pl.setGeometry(soilSpace) # creates soil space to stop roots from growing out of the soil

pl.initialize()
pl.simulate(simtime)

kz = 4.32e-2  # axial conductivity [cm3/day]
kr_root = 1.728e-4  # radial conductivity [1/day]
kr_stem = 0  # radial conductivity of stem  [1/day], set to 0
kr_leaf =  0.004 #  radial conductivity of leaf xylem membrane [1/day] different from stomatal conductance
p_top = -1000  # top soil total potential [cm]
p_bot = -900 # bot total potential [cm]

#envirenmental variables
RH = 0.5 # relative humidity [-]
TairC = 20 #air temperature [°C]
TairK = TairC + 273.15 #air temperature [K]
Q = 900e-6 #irradiance [mol quanta m-2 s-1]
cs = 350e-6 #co2 paartial pressure at leaf surface [mol mol-1]
es = 0.61078 * np.exp(17.27 * TairC / (TairC + 237.3))/10 # [hPa] #from FAO56
ea = es * RH
VPD = es - ea #vapor pressure deficit 
N_input = 4.4 #nitrogen content [%]


Rgaz=8.314 #J K-1 mol-1 = cm^3*MPa/K/mol
rho_h2o = 1 #g/cm3
Mh2o = 18.05 #g/mol
MPa2hPa = 10000
hPa2cm = 1.0197                     
p_a = np.log(RH)*Rgaz*rho_h2o*TairK/Mh2o *MPa2hPa*hPa2cm #air water potential [cm]


" Coupling to soil "
z_ = np.linspace(0,-35,36) # 0, -1, -2, ... -35
mids_ = 0.5*(z_[:-1] + z_[1:]) # 35
p_s = np.linspace(p_top, p_bot, 35)
soil_index = lambda x,y,z : max(int(np.floor(-z)),-1) #abovegroud nodes get index -1
pl.setSoilGrid(soil_index)
#print(p_s)
#raise exception
""" Set up organ conductivities"""

r = Leuning(pl) #Equivalent to XylemFluxPython() for whole plant
r.setKr([[kr_root],[kr_stem],[kr_leaf]]) #gmax will be changed by the leuning function 
r.setKx([[kz]])
r.airPressure = p_a 



#
# In the following 'cells = True' means total soil potentials 'p_s' is given for soil voxel
#
rx = r.solve_leuning( sim_time = simtime,sxx=p_s, cells = True, Qlight = Q,VPD = VPD,Tl = TairK,p_linit = p_s[0],
ci_init = cs*0.7, cs =cs, soil_k = [], N = N_input, log = False, verbose = True)


#fluxes1 = r.segFluxes(simtime, rx, p_s, cells = False)  # [cm3/day]

fluxes = r.radial_fluxes(simtime, rx, p_s,[], True)  # cm3/day
surfs = np.multiply(np.array(r.segLength()), 2*np.array(r.rs.radii)*np.pi)  # root segment side surface [cm2]
print("balances test ",sum(fluxes))
fluxes = np.divide(fluxes, surfs)  # we convert to [cm3/(cm2 day)]

""" plot results """
fig, ax = plt.subplots()
name = ["root", "stem", "leaf"]
color = ['tab:blue', 'tab:orange', 'tab:green']
for ndType in [2, 3, 4]:
    #nodes = self.get_nodes()
    y = r.get_nodes_organ_type(ndType)#coordinates
    x = np.array(rx)[r.get_nodes_index(ndType)]
    ax.scatter(x, y[:,2], c=color[ndType-2],  label=name[ndType-2],
               alpha=0.3, edgecolors='none')

ax.legend()
ax.grid(True)
plt.xlabel("Xylem pressure (cm)")
plt.ylabel("Depth (cm)")
plt.title("Xylem matric potential (cm)")
plt.show()


fig, ax = plt.subplots()
name = ["root", "stem", "leaf"]
color = ['tab:blue', 'tab:orange', 'tab:green']
for ndType in [2, 3, 4]:
    segIdx = r.get_segments_index(ndType)
    nodesy = segIdx + np.ones(segIdx.shape, dtype = np.int64)
    nodes = r.get_nodes()
    y = np.array(nodes)[nodesy]#coordinates
    x = fluxes[segIdx]
    ax.scatter(x, y[:,2], c=color[ndType-2],  label=name[ndType-2],
               alpha=0.3, edgecolors='none')

ax.legend()
ax.grid(True)
plt.xlabel("Fluxes (cm/day)")
plt.ylabel("Depth (cm)")
plt.title("water fluxes")
plt.show()


print("Transpiration", sum(x), "cm3/day")
print("Net assimilation", np.sum(r.An)*1e3, "mmol CO2 m-2 s-1")