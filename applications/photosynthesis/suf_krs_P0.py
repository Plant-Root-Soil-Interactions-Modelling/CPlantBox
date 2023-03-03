import sys;
sys.path.append("../..")#("/mnt/c/Users/mobil/CPlantBox"); 
sys.path.append("../../src/python_modules")#("/mnt/c/Users/mobil/CPlantBox/src/python_modules")
CPBdir = "../.."#"/mnt/c/Users/mobil/CPlantBox"

from xylem_flux import XylemFluxPython  # Python hybrid solver
import plantbox as pb
import vtk_plot as vp
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


""" Parameters """
# artificial shoot
# kr00 = np.array([[-154, 0.], [0, 1.e-12], [1e20, 1.e-12]])
# kz00 = np.array([[-154, 0.], [0, 1.], [1e20, 1.]])

# kr_r0 = list(pd.read_csv('/mnt/c/Users/mobil/CPlantBox_test_files/kr/kr_Tap.csv').k)
# kr_age_r0 = list(pd.read_csv('/mnt/c/Users/mobil/CPlantBox_test_files/kr/kr_Tap.csv').age)
# kr_r1 = list(pd.read_csv('/mnt/c/Users/mobil/CPlantBox_test_files/kr/kr_Basal.csv').k)
# kr_age_r1 = list(pd.read_csv('/mnt/c/Users/mobil/CPlantBox_test_files/kr/kr_Basal.csv').age)
# kr_r2 = list(pd.read_csv('/mnt/c/Users/mobil/CPlantBox_test_files/kr/kr_SBR.csv').k)
# kr_age_r2 = list(pd.read_csv('/mnt/c/Users/mobil/CPlantBox_test_files/kr/kr_SBR.csv').age)
# kr_r3 = list(pd.read_csv('/mnt/c/Users/mobil/CPlantBox_test_files/kr/kr_LAT_A.csv').k)
# kr_age_r3 = list(pd.read_csv('/mnt/c/Users/mobil/CPlantBox_test_files/kr/kr_LAT_A.csv').age)
# kr_r4 = list(pd.read_csv('/mnt/c/Users/mobil/CPlantBox_test_files/kr/kr_LAT.csv').k)
# kr_age_r4 = list(pd.read_csv('/mnt/c/Users/mobil/CPlantBox_test_files/kr/kr_LAT.csv').age)

# kz_r0 = list(pd.read_csv('/mnt/c/Users/mobil/CPlantBox_test_files/kx/kx_Tap.csv').k)
# kz_age_r0 = list(pd.read_csv('/mnt/c/Users/mobil/CPlantBox_test_files/kx/kx_Tap.csv').age)
# kz_r1 = list(pd.read_csv('/mnt/c/Users/mobil/CPlantBox_test_files/kx/kx_Basal.csv').k)
# kz_age_r1 = list(pd.read_csv('/mnt/c/Users/mobil/CPlantBox_test_files/kx/kx_Basal.csv').age)
# kz_r2 = list(pd.read_csv('/mnt/c/Users/mobil/CPlantBox_test_files/kx/kx_SBR.csv').k)
# kz_age_r2 = list(pd.read_csv('/mnt/c/Users/mobil/CPlantBox_test_files/kx/kx_SBR.csv').age)
# kz_r3 = list(pd.read_csv('/mnt/c/Users/mobil/CPlantBox_test_files/kx/kx_LAT_A.csv').k)
# kz_age_r3 = list(pd.read_csv('/mnt/c/Users/mobil/CPlantBox_test_files/kx/kx_LAT_A.csv').age)
# kz_r4 = list(pd.read_csv('/mnt/c/Users/mobil/CPlantBox_test_files/kx/kx_LAT.csv').k)
# kz_age_r4 = list(pd.read_csv('/mnt/c/Users/mobil/CPlantBox_test_files/kx/kx_LAT.csv').age)


np.set_printoptions(threshold=sys.maxsize)
kz = 4.32e-2  # [cm^3/day]
kr = 1.728e-4  # [1/day]

simtime = 28  # [day] for task b

""" root system """
rs = pb.MappedRootSystem()
rs.setGeometry(pb.SDF_PlantBox(1.e6, 1.e6, 1.e6))  # not allowed to grow upwards out of soil
p_s = np.linspace(-500, -200, 3001)  #  -200.*np.ones((2001, 1))   # 3 meter down, from -200 to -500, resolution in mm
soil_index = lambda x, y, z: int(-10 * z)  # maps to p_s (hydrostatic equilibirum)
rs.setSoilGrid(soil_index)

path = "../../modelparameter/plant/"#"/mnt/c/Users/mobil/CPlantBox_test_files/params/"
# path = '/mnt/c/Users/mobil/CPlantBox/modelparameter/rootsystem/'
name = "P0_rs"  # "Glycine_max"  # "Anagallis_femina_Leitner_2010"  # Zea_mays_1_Leitner_2010
rs.readParameters(path + name + ".xml")
#rs.getRootSystemParameter().seedPos.z = -3
rs.setSeed(2)  # random
print(path + name + ".xml")

#print("non ",rs.getNumberOfNodes())
#raise Exception
r = XylemFluxPython(rs)#getNumberOfNodes
r.rs.initialize(True)
print(r.rs.getNumberOfNewNodes(),r.rs.getNumberOfNodes(),len(r.get_nodes()))
print(r.rs.getNumberOfOrgans())
#raise Exception
print(r.rs.getNumberOfNewNodes())
print(r.get_segments())
print(r.get_subtypes())
#raise Exception
r.rs.simulate(simtime, True)
print(r.rs.getNumberOfNewNodes(),r.rs.getNumberOfNodes(),len(r.get_nodes()))
#raise Exception
""" set up xylem parameters """

r.setKr([[kr]])#[kr00[:, 1],kr_r0, kr_r1, kr_r2, kr_r2, kr_r3, kr_r4], [kr00[:, 0],kr_age_r0,kr_age_r1,kr_age_r2,kr_age_r2,kr_age_r3,kr_age_r4])
r.setKx([[kz]])#[kz00[:, 1],kz_r0,kz_r1,kz_r2,kz_r2,kz_r3,kz_r4], [kz00[:, 0],kz_age_r0,kz_age_r1,kz_age_r2,kz_age_r2,kz_age_r3,kz_age_r4])

""" for debugging """
print(r.get_subtypes())
r.test()
r.plot_conductivities()
shoot_segs = rs.getShootSegments()
print("Shoot segments", [str(s) for s in shoot_segs])
print("Shoot type", rs.subTypes[0])
#print(r.rs.organTypes,len(r.rs.organTypes),len(r.rs.subTypes))
raise Exception("finished")

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

