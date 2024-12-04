import sys;

CPBdir = "../.."
sys.path.append(CPBdir + "/src");
sys.path.append(CPBdir);
sys.path.append("../../..");
sys.path.append("..");
sys.path.append(CPBdir + "/src/python_modules");
sys.path.append("../build-cmake/cpp/python_binding/")  # dumux python binding
sys.path.append("../../build-cmake/cpp/python_binding/")
sys.path.append("../modules/")  # python wrappers

import helpuqrMasterCopy1
import importlib
import pandas as pd

importlib.reload(helpuqrMasterCopy1)
from helpuqrMasterCopy1 import *
from helpuqrMasterCopy1 import setKrKx_xylem
# reload(helpuqrMasterCopy1)
import numpy as np


def weather(simDuration, hp, condition):
    vgSoil = [0.059, 0.45, 0.00644, 1.503, 1]
    loam = [0.08, 0.43, 0.04, 1.6, 50]
    Qnigh = 0;
    Qday = 960e-6  # 458*2.1
    if (condition == "wet"):
        Tnigh = 15.8;
        Tday = 22
        # Tnigh = 13; Tday = 20.7
        # specificHumidity = 0.0097
        RHday = 0.60;
        RHnigh = 0.88
        Pair = 1010.00  # hPa
        thetaInit = 40 / 100  # 15.59/100#
        cs = 350e-6
    elif condition == "dry":
        Tnigh = 20.7;
        Tday = 30.27
        # Tnigh = 15.34; Tday = 23.31
        # specificHumidity = 0.0097# 0.0111
        RHday = 0.3;
        RHnigh = 0.78
        Pair = 1070.00  # hPa
        thetaInit = 28 / 100  # 10.47/100#
        cs = 350e-6
    else:
        print("condition", condition)
        raise Exception("condition not recognised")
    coefhours = 1# sinusoidal(simDuration)
    RH_ = RHnigh + (RHday - RHnigh) * coefhours
    TairC_ = Tnigh + (Tday - Tnigh) * coefhours
    Q_ = Qnigh + (Qday - Qnigh) * coefhours
    # co2 paartial pressure at leaf surface (mol mol-1)
    # 390, 1231
    # RH = 0.5 # relative humidity
    es = 6.112 * np.exp((17.67 * TairC_) / (TairC_ + 243.5))
    ea = es * RH_  # qair2ea(specificHumidity,  Pair)
    assert ea < es
    # RH = ea/es
    assert ((RH_ > 0) and (RH_ < 1))
    bl_thickness = 1 / 1000  # m
    diffusivity = 2.5e-5  # m2/sfor 25°C
    rbl = bl_thickness / diffusivity  # s/m 13
    # cs = 350e-6
    Kcanopymean = 1e-1  # m2/s
    meanCanopyL = (2 / 3) * hp / 2
    rcanopy = meanCanopyL / Kcanopymean
    windSpeed = 2  # m/s
    zmzh = 2  # m
    karman = 0.41  # [-]

    rair = 1
    if hp > 0:
        rair = np.log((zmzh - (2 / 3) * hp) / (0.123 * hp)) * np.log(
            (zmzh - (2 / 3) * hp) / (0.1 * hp)) / (karman * karman * windSpeed)
        # print()
        # raise Exception

    pmean = theta2H(vgSoil, thetaInit)

    weatherVar = {'TairC': TairC_, 'TairK': TairC_ + 273.15, 'Pair': Pair, "es": es,
                  'Qlight': Q_, 'rbl': rbl, 'rcanopy': rcanopy, 'rair': rair, "ea": ea,
                  'cs': cs, 'RH': RH_, 'p_mean': pmean, 'vg': loam}
    # print("Env variables at", round(simDuration//1),"d",round((simDuration%1)*24),"hrs :\n", weatherVar)
    return weatherVar


def resistance2conductance(resistance, r, weatherX):
    resistance = resistance * (1 / 100)  # [s/m] * [m/cm] = [s/cm]
    resistance = resistance * r.R_ph * weatherX[
        "TairK"] / r.Patm  # [s/cm] * [K] * [hPa cm3 K−1 mmol−1] * [hPa] = [s] * [cm2 mmol−1]
    resistance = resistance * (1000) * (
                1 / 10000)  # [s cm2 mmol−1] * [mmol/mol] * [m2/cm2] = [s m2 mol−1]
    return 1 / resistance


def initPlant(simDuration, condition):
    weatherInit = weather(0, 0, condition)
    # simDuration = 25 # [day] init simtime
    # spellDuration = 5
    simMax = 26  # simStartSim+ spellDuration
    depth = 40
    dt = 1 / 24  # 10min
    verbose = True

    # plant system
    pl = pb.MappedPlant(seednum=2)  # pb.MappedRootSystem() #pb.MappedPlant()
    # pl2 = pb.MappedPlant(seednum = 2) #pb.MappedRootSystem() #pb.MappedPlant()
    path = CPBdir + "/modelparameter/structural/plant/"
    name = "P0_plant_new"  # "Triticum_aestivum_adapted_2021"#

    pl.readParameters( name + ".xml")

    # raise Exception
    sdf = pb.SDF_PlantBox(np.inf, np.inf, depth)

    pl.setGeometry(sdf)  # creates soil space to stop roots from growing out of the soil

    pl.initialize(verbose=True)  # , stochastic = False)
    
    parameters = pl.getOrganRandomParameter(pb.stem)
    for p in parameters:
        print(p)
    parameters = pl.getOrganRandomParameter(pb.root)
    for p in parameters:
        print(p)
    parameters = pl.getOrganRandomParameter(pb.leaf)
    for p in parameters:
        print(p)
    
    pl.simulate(simDuration, False)  # , "outputpm15.txt")

    # raise Exception
    """ Coupling to soil """

    min_b = [-3. / 2, -12. / 2, -41.]  # distance between wheat plants
    max_b = [3. / 2, 12. / 2, 0.]
    rez = 0.5
    cell_number = [int(6 * rez), int(24 * rez), int(40 * rez)]  # 1cm3?
    layers = depth;
    soilvolume = (depth / layers) * 3 * 12
    k_soil = []
    initial = weatherInit["p_mean"]  # mean matric potential [cm] pressure head

    p_mean = initial
    p_bot = p_mean + depth / 2
    p_top = initial - depth / 2
    sx = np.linspace(p_top, p_bot, depth)
    picker = lambda x, y, z: max(int(np.floor(-z)), -1)
    sx_static_bu = sx
    pl.setSoilGrid(picker)  # maps segment

    """ Parameters phloem and photosynthesis """
    r = PhloemFluxPython(pl, psiXylInit=min(sx),
                         ciInit=weatherInit["cs"] * 0.5)  # XylemFluxPython(pl)#
    # r2 = PhloemFluxPython(#pl2,psiXylInit = min(sx),ciInit = weatherInit["cs"]*0.5) #XylemFluxPython(pl)#

    r = setKrKx_phloem(r)

    r.oldciEq = True

    
    r.g0 = 8e-3
    r.VcmaxrefChl1 = 1.1#1.28
    r.VcmaxrefChl2 = 4#8.33
    r.a1 = 0.6 / 0.4  # 0.7/0.3#0.6/0.4 #ci/(cs - ci) for ci = 0.6*cs
    r.a3 = 2
    r.alpha =0.4# 0.2#/2
    r.theta = 0.6#0.9#/2
    r.k_meso = 1e-3  # 1e-4
    r.setKrm2([[2e-5]])
    r.setKrm1([[10e-2]])  # ([[2.5e-2]])
    r.setRhoSucrose([[0.51], [0.65], [0.56]])  # 0.51
    # ([[14.4,9.0,0,14.4],[5.,5.],[15.]])
    rootFact = 2
    r.setRmax_st(
        [[2.4 * rootFact, 1.5 * rootFact, 0.6 * rootFact, 2.4 * rootFact], [2., 2.],
         [8.]])  # 6.0#*6 for roots, *1 for stem, *24/14*1.5 for leaves
    # r.setRmax_st([[12,9.0,6.0,12],[5.,5.],[15.]])
    r.KMrm = 0.1  # VERY IMPORTANT TO KEEP IT HIGH
    # r.exud_k = np.array([2.4e-4])#*10#*(1e-1)
    # r.k_gr = 1#0
    r.sameVolume_meso_st = False
    r.sameVolume_meso_seg = True
    r.withInitVal = True
    r.initValST = 0.  # 0.6#0.0
    r.initValMeso = 0.  # 0.9#0.0
    r.beta_loading = 0.6
    r.Vmaxloading = 0.05  # mmol/d, needed mean loading rate:  0.3788921068507634
    r.Mloading = 0.2
    r.Gr_Y = 0.8
    r.CSTimin = 0.4
    r.surfMeso = 0.0025
    r.leafGrowthZone = 2  # cm
    r.StemGrowthPerPhytomer = True  #
    r.psi_osmo_proto = -10000 * 1.0197  # schopfer2006
    r.fwr = 0

    r.cs = weatherInit["cs"]

    # r.r_forPhloem(24/14*1.5, 4)
    # r.r_forPhloem(24/14, 3)
    # r.r_forPhloem(6, 2) #because roots will have high C_ST less often
    r.expression = 6
    r.update_viscosity = True
    r.solver = 1
    r.atol = 1e-10
    r.rtol = 1e-6
    # r.doNewtonRaphson = False;r.doOldEq = False
    SPAD = 31.0
    chl_ = (0.114 * (SPAD ** 2) + 7.39 * SPAD + 10.6) / 10
    r.Chl = np.array([chl_])
    r.Csoil = 1e-4

    hp = max([tempnode[2] for tempnode in r.get_nodes()]) / 100

    weatherX = weather(simDuration, hp, condition)
    r.Patm = weatherX["Pair"]
    ##resistances
    r.g_bl = resistance2conductance(weatherX["rbl"], r, weatherX) / r.a2_bl
    r.g_canopy = resistance2conductance(weatherX["rcanopy"], r, weatherX) / r.a2_canopy
    r.g_air = resistance2conductance(weatherX["rair"], r, weatherX) / r.a2_air
    r.sh = 4e-4
    r.fwr = 0  # 0.001
    r.shmesophyll = 4e-4
    r.fwrmesophyll = 0  # 0.001
    r.p_lcritmesophyll = -100000
    # r.gm = 0.05
    r.g0 = 8e-6
    r.alternativeAn = False
    r.limMaxErr = 1 / 100;

    r.Qlight = weatherX["Qlight"]  # ; TairC = weatherX["TairC"] ; text = "night"
    kr_l = 3.83e-5
    r = setKrKx_xylem(weatherX["TairC"], weatherX["RH"], r, kr_l)
    kr_l = 3.83e-5
    # r.setKr_meso([kr_l])
    r.es = weatherX["es"]
    return r, weatherX, sx


Q = np.arange(1e-6, 1000e-6, 100e-6)  # mol quanta m-2 s-1 light, example from leuning1995
TairC = np.arange(-10, 55, 5)
cs = np.arange(0e-6, 500e-6, 50e-6)  # mol mol-1

simDuration = 10
rinit = initPlant(simDuration, condition="wet")
r = rinit[0]
weatherX = rinit[1]
rinit = initPlant(simDuration, condition="wet")
rC4 = rinit[0]
print("r.PhotoType, pb.C3, pb.C4",r.PhotoType,rC4.PhotoType)
rC4.PhotoType=pb.C4


sx = rinit[2]  # np.linspace(-100, -20000, 100)
p_errors = []
dt = 1 / 24
kr_l = 3.83e-5
r = setKrKx_xylem(weatherX["TairC"], weatherX["RH"], r, kr_l)
rC4 = setKrKx_xylem(weatherX["TairC"], weatherX["RH"], rC4, kr_l)
#weatherX = weather(simDuration, hp, "wet")
r.Patm = weatherX["Pair"]
rC4.Patm = weatherX["Pair"]
directoryN = "/testPhoto2/"
variables = [Q,TairC,  cs]
#variables = [[],[],  cs]
resultsAnC3 = [[], [], []]
resultsAnC4 = [[], [], []]
print("to do",[len(a) for a in variables])
for i in range(len(variables)):
    An = []
    AnC4 = []
    Q_input = 900e-6
    cs_input = 350e-6
    Tair_input = 25
    var = [Q_input, Tair_input, cs_input]
    for j in range(len(variables[i])):
        print(i, j)
        var[i] = variables[i][j]
        
        weatherX["Qlight"] = var[0]
        weatherX["TairC"] = var[1]
        weatherX["TairK"] = weatherX["TairC"] + 273.15
        weatherX["cs"] = var[2]
        r.cs = weatherX["cs"]
        rC4.cs = weatherX["cs"]
        hp = max([tempnode[2] for tempnode in r.get_nodes()]) / 100  # maxnode canopy [m]
        # print([tempnode[2] for tempnode in r.get_nodes()], hp)

        
        ##resistances
        r.g_bl = resistance2conductance(weatherX["rbl"], r, weatherX) / r.a2_bl
        r.g_canopy = resistance2conductance(weatherX["rcanopy"], r, weatherX) / r.a2_canopy
        r.g_air = resistance2conductance(weatherX["rair"], r, weatherX) / r.a2_air
        rC4.g_bl = resistance2conductance(weatherX["rbl"], rC4, weatherX) / rC4.a2_bl
        rC4.g_canopy = resistance2conductance(weatherX["rcanopy"], rC4, weatherX) / rC4.a2_canopy
        rC4.g_air = resistance2conductance(weatherX["rair"], rC4, weatherX) / rC4.a2_air

        r.Qlight = weatherX["Qlight"]  # ; TairC = weatherX["TairC"] ; text = "night"
        rC4.Qlight = weatherX["Qlight"]  
        
        
        
        #
        # kr_l  = 3.83e-7
        # r.setKr_meso([kr_l])
        print("r.PhotoType, pb.C3, pb.C4",r.PhotoType,rC4.PhotoType)
        r.solve_photosynthesis(sim_time_=simDuration, sxx_=sx, cells_=True,
                               ea_=weatherX["ea"], es_=weatherX["es"],
                               verbose_=False, doLog_=False, TairC_=weatherX["TairC"],
                               outputDir_="./results" + directoryN)
        rC4.solve_photosynthesis(sim_time_=simDuration, sxx_=sx, cells_=True,
                               ea_=weatherX["ea"], es_=weatherX["es"],
                               verbose_=False, doLog_=False, TairC_=weatherX["TairC"],
                               outputDir_="./results" + directoryN)
        idsC4 = np.where(np.array(rC4.ci)> 0)
        idsC3 = np.where(np.array(r.ci)> 0)
        assert (idsC4[0] == idsC3[0]).all()
        print('Qlight:{}, T:{}, cs:{},An:{}, Vj:{}, Vp:{},Vc:{}, Rd:{}, Rd_ref:{}'.format(var[0]*1e6, var[1], var[2]*1e6,
                np.mean(np.array(rC4.An)[idsC3])*1e6,np.mean(np.array(rC4.Vj)[idsC3])*1e6,
                np.mean(np.array(rC4.Vp)[idsC3])*1e6,np.mean(np.array(rC4.Vc)[idsC3])*1e6,
                np.mean(np.array(rC4.Rd)[idsC3])*1e6,np.mean(np.array(rC4.Rd_ref)[idsC3])*1e6)) 
        print('An:{}, Vj:{}, Vp:{},Vc:{}, Rd:{}, Rd_ref:{}'.format(
                np.mean(np.array(r.An)[idsC3])*1e6,np.mean(np.array(r.Vj)[idsC3])*1e6,
                np.mean(np.array(r.Vp)[idsC3])*1e6,np.mean(np.array(r.Vc)[idsC3])*1e6,
                np.mean(np.array(r.Rd)[idsC3])*1e6,np.mean(np.array(r.Rd_ref)[idsC3])*1e6))
        #print(max(rC4.deltagco2),max(rC4.ci),max(rC4.gco2)) 
        #print(var,np.mean(rC4.Vc),np.mean(rC4.Vj),np.mean(rC4.Vp),np.mean(rC4.Rd),np.mean(rC4.Rd_ref))                       
        An.append(np.mean(np.array(r.An)[idsC3]) * 1e6)
        AnC4.append(np.mean(np.array(rC4.An)[idsC3]) * 1e6)
    resultsAnC3[i] = An
    resultsAnC4[i] = AnC4

resultsAnC3 = pd.DataFrame(resultsAnC3)
resultsAnC4 = pd.DataFrame(resultsAnC4)
variables = pd.DataFrame(variables)
resultsAnC3.to_csv("resultsAnC3.csv", index=False)
resultsAnC4.to_csv("resultsAnC4.csv", index=False)
variables.to_csv("variables.csv", index=False)

isInSoil = np.array([rC4.plant.seg2cell[si] >= 0 for si in range(len(rC4.plant.segments))])
ana = pb.SegmentAnalyser(rC4.plant.mappedSegments())
ana.addData("isInSoil", isInSoil)
ana.write("testphoto3.vtp", ["organType", "subType", "isInSoil"]) 
