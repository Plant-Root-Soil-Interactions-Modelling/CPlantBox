""" water movement within the root (static soil) """
directoryN = "/21to28wet/"
import sys; 
CPBdir = "../../../../CPB2506/CPlantBox"
sys.path.append(CPBdir+"/src");
sys.path.append(CPBdir);
sys.path.append("../../..");sys.path.append(".."); 
sys.path.append(CPBdir+"/src/python_modules");
sys.path.append("../build-cmake/cpp/python_binding/") # dumux python binding
sys.path.append("../../build-cmake/cpp/python_binding/")
sys.path.append("../modules/") # python wrappers 

#import matplotlib
#matplotlib.use('AGG') #otherwise "import matplotlib.pyplot" hangs

from rosi_richards import RichardsSP  # C++ part (Dumux binding)
from richards import RichardsWrapper  # Python part
from phloem_flux import PhloemFluxPython  # Python hybrid solver
#from Leuning_speedup import Leuning #about 0.7 for both
#from photosynthesis_cython import Leuning
import plantbox as pb
#from plantbox import Photosynthesis as ph
#import vtk_plot as vp
import math
import os
import numpy as np
import vtk_plot as vp
#import matplotlib.pyplot as plt
from datetime import datetime, timedelta


def delete_verbose():
    home_dir = os.getcwd()
    dir_name = ""#"/results"
    dir_name2 = home_dir + dir_name
    test = os.listdir(dir_name2)
    for item in test:
        if item.startswith("outphoto"):
            os.remove(os.path.join(dir_name2, item))
        if item.startswith("fluxesphoto"):
            os.remove(os.path.join(dir_name2, item))
        if item.startswith("loopphoto"):
            os.remove(os.path.join(dir_name2, item))
        if item.startswith("errphoto"):
            os.remove(os.path.join(dir_name2, item))
delete_verbose()
#       qr, qs, alpha, n, ks (in [cm/d])
#vgSoil = [0.059, 0.45, 0.00644, 1.503, 1]
def theta2H(vg,theta):#(-) to cm
    thetar =vg[0]# 0.059
    thetas = vg[1]#0.445
    alpha = vg[2]#0.00644
    n = vg[3]#1.503
    nrev = 1/(1-1/n)
    H =-(((( (thetas - thetar)/(theta - thetar))**nrev) - 1)**(1/n))/alpha
    return(H)#cm

def sinusoidal(t):
    return (np.sin(np.pi*t*2)+1)/2

def qair2rh(qair, es_, press):
    e =qair * press / (0.378 * qair + 0.622)
    rh = e / es_
    rh=max(min(rh, 1.),0.)
    return rh


def weather(simDuration):
    vgSoil = [0.059, 0.45, 0.00644, 1.503, 1]
    
    Qmin = 0; Qmax = 960e-6 #458*2.1
    Tmin = 15.8; Tmax = 22
    specificHumidity = 0.0097
    Pair = 1010.00 #hPa
    thetaInit = 30/100

    coefhours = sinusoidal(simDuration)
    TairC_ = Tmin + (Tmax - Tmin) * coefhours
    Q_ = Qmin + (Qmax - Qmin) * coefhours
    cs = 850e-6 #co2 paartial pressure at leaf surface (mol mol-1)
    #RH = 0.5 # relative humidity
    es =  6.112 * np.exp((17.67 * TairC_)/(TairC_ + 243.5))
    RH = qair2rh(specificHumidity, es, Pair)
    
    pmean = theta2H(vgSoil, thetaInit)
    
    weatherVar = {'TairC' : TairC_,
                    'Qlight': Q_,
                    'cs':cs, 'RH':RH, 'p_mean':pmean, 'vg':vgSoil}
    print("Env variables at", round(simDuration//1),"d",round((simDuration%1)*24),"hrs :\n", weatherVar)
    return weatherVar

def div0(a, b, c):        
    return np.divide(a, b, out=np.full(len(a), c), where=b!=0)
    
def div0f(a, b, c):    
    if b != 0:
        return a/b
    else:
        return a/c
        
def write_file_array(name, data):
    name2 = 'results'+ directoryN+ name+ '_15pm.txt'
    with open(name2, 'a') as log:
        log.write(','.join([num for num in map(str, data)])  +'\n')

def write_file_float(name, data):
    name2 = 'results' + directoryN+  name+ '_15pm.txt'
    with open(name2, 'a') as log:
        log.write(repr( data)  +'\n')

home_dir = os.getcwd()
dir_name = "/results"+directoryN
dir_name2 = home_dir + dir_name
test = os.listdir(dir_name2)
for item in test:
    if item.endswith("_15pm.txt"):
        os.remove(os.path.join(dir_name2, item))
    if item.endswith("uqr15.vtk"):
        os.remove(os.path.join(dir_name2, item))
    if item.startswith("15pm"):
        os.remove(os.path.join(dir_name2, item))

def setKrKx_xylem(TairC, RH): #inC
    #mg/cm3
    hPa2cm = 1.0197
    dEauPure = (999.83952 + TairC * (16.952577 + TairC * 
        (- 0.0079905127 + TairC * (- 0.000046241757 + TairC * 
        (0.00000010584601 + TairC * (- 0.00000000028103006)))))) /  (1 + 0.016887236 * TairC)
    siPhi = (30 - TairC) / (91 + TairC)
    siEnne=0
    mu =  pow(10, (- 0.114 + (siPhi * (1.1 + 43.1 * pow(siEnne, 1.25) )))) 
    mu = mu /(24*60*60)/100/1000; #//mPa s to hPa d, 1.11837e-10 hPa d for pure water at 293.15K
    mu = mu * hPa2cm #hPa d to cmh2o d 

    #number of vascular bundles
    VascBundle_leaf = 32
    VascBundle_stem = 52
    VascBundle_root = 1 #valid for all root type
            
    #radius of xylem type^4 * number per bundle
    rad_x_l_1   = (0.0015 **4) * 2; rad_x_l_2   = (0.0005 **4) * 2   
    rad_x_s_1   = (0.0017 **4) * 3; rad_x_s_2   = (0.0008 **4) * 1     
    rad_x_r0_1  = (0.0015 **4) * 4    
    rad_x_r12_1 = (0.00041**4) * 4; rad_x_r12_2 = (0.00087**4) * 1
    rad_x_r3_1  = (0.00068**4) * 1      

    # axial conductivity [cm^3/day]        
    kz_l  = VascBundle_leaf *(rad_x_l_1 + rad_x_l_2)    *np.pi /(mu * 8)  
    kz_s  = VascBundle_stem *(rad_x_s_1 + rad_x_s_2)    *np.pi /(mu * 8) 
    kz_r0 = VascBundle_root * rad_x_r0_1                *np.pi /(mu * 8)  
    kz_r1 = VascBundle_root *(rad_x_r12_1 + rad_x_r12_2)*np.pi /(mu * 8) 
    kz_r2 = VascBundle_root *(rad_x_r12_1 + rad_x_r12_2)*np.pi /(mu * 8)  
    kz_r3 = VascBundle_root * rad_x_r3_1                *np.pi /(mu * 8) # 4.32e-1

    #radial conductivity [1/day],
    kr_l  = 3.83e-4 * hPa2cm# init: 3.83e-4 cm/d/hPa
    kr_s  = 0.#1.e-20  * hPa2cm # set to almost 0
    kr_r0 = 6.37e-5 * hPa2cm 
    kr_r1 = 7.9e-5  * hPa2cm 
    kr_r2 = 7.9e-5  * hPa2cm  
    kr_r3 = 6.8e-5  * hPa2cm 
    l_kr = 0.8 #cm
    r.setKr([[kr_r0,kr_r1,kr_r2,kr_r0],[kr_s,kr_s ],[kr_l]], kr_length_=l_kr) 
    r.setKx([[kz_r0,kz_r1,kz_r2,kz_r0],[kz_s,kz_s ],[kz_l]])
    
    
    Rgaz=8.314 #J K-1 mol-1 = cm^3*MPa/K/mol
    rho_h2o = dEauPure/1000#g/cm3
    Mh2o = 18.05 #g/mol
    MPa2hPa = 10000
    hPa2cm = 1/0.9806806
    #log(-) * (cm^3*MPa/K/mol) * (K) *(g/cm3)/ (g/mol) * (hPa/MPa) * (cm/hPa) =  cm                      
    p_a = np.log(RH) * Rgaz * rho_h2o * (TairC + 273.15)/Mh2o * MPa2hPa * hPa2cm

    r.psi_air = p_a #*MPa2hPa #used only with xylem

    
def setKrKx_phloem(): #inC

    #number of vascular bundles
    VascBundle_leaf = 32
    VascBundle_stem = 52
    VascBundle_root = 1 #valid for all root type
            
    #numPerBundle
    numL = 18
    numS = 21
    numr0 = 33
    numr1 = 25
    numr2 = 25
    numr3 = 1
        
    #radius of phloem type^4 * number per bundle
    rad_s_l   = numL* (0.00025 **4)# * 2; rad_x_l_2   = (0.0005 **4) * 2   
    rad_s_s   = numS *(0.00019 **4) #* 3; rad_x_s_2   = (0.0008 **4) * 1     
    rad_s_r0  = numr0 *(0.00039 **4) #* 4    
    rad_s_r12 = numr1*(0.00035**4) #* 4; rad_x_r12_2 = (0.00087**4) * 1
    rad_s_r3  = numr3 *(0.00068**4) #* 1      

    # axial conductivity [cm^3/day] , mu is added later as it evolves with CST  
    beta = 0.9 #Thompson 2003a
    kz_l   = VascBundle_leaf * rad_s_l   * np.pi /8 * beta  
    kz_s   = VascBundle_stem * rad_s_s   * np.pi /8 * beta
    kz_r0  = VascBundle_root * rad_s_r0  * np.pi /8 * beta
    kz_r12 = VascBundle_root * rad_s_r12 * np.pi /8 * beta
    kz_r3  = VascBundle_root * rad_s_r3  * np.pi /8 * beta
    
    #print([[kz_r0,kz_r12,kz_r12,kz_r3],[kz_s,kz_s ],[kz_l]])
    #raise Exception
    #radial conductivity [1/day],
    kr_l  = 0.#3.83e-4 * hPa2cm# init: 3.83e-4 cm/d/hPa
    kr_s  = 0.#1.e-20  * hPa2cm # set to almost 0
    kr_r0 = 5e-2
    kr_r1 = 5e-2
    kr_r2 = 5e-2
    kr_r3 = 5e-2
    l_kr = 0.8 #cm
    
    r.setKr_st([[kr_r0,kr_r1 ,kr_r2 ,kr_r0],[kr_s,kr_s ],[kr_l]] , kr_length_= l_kr)
    r.setKx_st([[kz_r0,kz_r12,kz_r12,kz_r0],[kz_s,kz_s ],[kz_l]])
    
    a_ST = [[0.00039,0.00035,0.00035,0.00039 ],[0.00019,0.00019],[0.00025]]
    Across_s_l   = numL*VascBundle_leaf *(a_ST[2][0]**2)*np.pi# (0.00025 **2)# * 2; rad_x_l_2   = (0.0005 **4) * 2   
    Across_s_s   = numS *VascBundle_stem * (a_ST[1][0]**2)*np.pi#(0.00019 **2) #* 3; rad_x_s_2   = (0.0008 **4) * 1     
    Across_s_r0  = numr0 *VascBundle_root * (a_ST[0][0]**2)*np.pi#(0.00039 **2) #* 4    
    Across_s_r12 = numr1*VascBundle_root * (a_ST[0][1]**2)*np.pi#(0.00035**2) #* 4; rad_x_r12_2 = (0.00087**4) * 1
    Across_s_r3  =  numr3 *VascBundle_root *(a_ST[0][2]**2)*np.pi# (0.00068**2) #* 1    
    
    Perimeter_s_l   = numL*VascBundle_leaf *(a_ST[2][0])* 2 * np.pi# (0.00025 **2)# * 2; rad_x_l_2   = (0.0005 **4) * 2   
    Perimeter_s_s   = numS *VascBundle_stem * (a_ST[1][0])* 2 * np.pi#(0.00019 **2) #* 3; rad_x_s_2   = (0.0008 **4) * 1     
    Perimeter_s_r0  = numr0 *VascBundle_root * (a_ST[0][0])* 2 * np.pi#(0.00039 **2) #* 4    
    Perimeter_s_r12 = numr1*VascBundle_root * (a_ST[0][1])* 2 * np.pi#(0.00035**2) #* 4; rad_x_r12_2 = (0.00087**4) * 1
    Perimeter_s_r3  =  numr3 *VascBundle_root *(a_ST[0][2])* 2 * np.pi# (0.00068**2) #* 1  
    #print(a_ST[2][0],a_ST[1][0],a_ST[0][0],a_ST[0][1],a_ST[0][2])
    #r.a_ST = a_ST #to check for water equilibrium assumption
    #tot surface/np.pi of sieve tube  (np.pi added after)
    #r.a_ST_eqs = [[rad_s_r0,rad_s_r12,rad_s_r12,rad_s_r0],[rad_s_s,rad_s_s],[rad_s_l]]
    r.setAcross_st([[Across_s_r0,Across_s_r12,Across_s_r12,Across_s_r0],[Across_s_s,Across_s_s],[Across_s_l]])
    #actually, don t use perimeter currently
    #r.setPerimeter_st([[0,Perimeter_s_r0,Perimeter_s_r12,Perimeter_s_r12,Perimeter_s_r0],[0,Perimeter_s_s,Perimeter_s_s],[0,Perimeter_s_l]])
    
""" Parameters """

weatherInit = weather(0)
simInit = 21
simDuration = simInit # [day] init simtime
simMax =simInit + 7
depth = 60
dt = 1/24 #1h
verbose = True

# plant system 
pl = pb.MappedPlant(seednum = 2) #pb.MappedRootSystem() #pb.MappedPlant()
pl2 = pb.MappedPlant(seednum = 2) #pb.MappedRootSystem() #pb.MappedPlant()
path = CPBdir+"/modelparameter/plant/"
name = "Triticum_aestivum_adapted_2021"#"wheat_uqr15" #"manyleaves"## "Anagallis_femina_Leitner_2010"  # Zea_mays_1_Leitner_2010

pl.readParameters(path + name + ".xml")
pl2.readParameters(path + name + ".xml")


#raise Exception
sdf = pb.SDF_PlantBox(np.Inf, np.Inf, depth )

pl.setGeometry(sdf) # creates soil space to stop roots from growing out of the soil
pl2.setGeometry(sdf) # creates soil space to stop roots from growing out of the soil


pl.initialize(verbose = True)#, stochastic = False)
pl.simulate(simDuration, False)#, "outputpm15.txt")
pl2.initialize(verbose = True)#, stochastic = False)
pl2.simulate(simDuration, False)#, "outputpm15.txt")
ot_ =np.array(pl.organTypes)
segments_ = np.array(pl.segLength())
print(len(segments_), sum(segments_))
print(len(segments_[np.where(ot_ ==2)]), sum(segments_[np.where(ot_ ==2)]))
print(len(segments_[np.where(ot_ ==3)]), sum(segments_[np.where(ot_ ==3)]))
print(len(segments_[np.where(ot_ ==4)]), sum(segments_[np.where(ot_ ==4)]))
print(np.array(pl.leafBladeSurface)[np.where(ot_ ==4)],
    len(pl.leafBladeSurface))
#raise Exception
""" Coupling to soil """



min_b = [-3./2, -12./2, -61.]#distance between wheat plants
max_b = [3./2, 12./2, 0.]
cell_number = [6, 24, 61]#1cm3? 
layers = depth; soilvolume = (depth / layers) * 3 * 12
k_soil = []
initial = weatherInit["p_mean"]#mean matric potential [cm] pressure head

s = RichardsWrapper(RichardsSP())
s.initialize()
periodic = True
s.createGrid(min_b, max_b, cell_number, periodic)  # [cm]
s.setHomogeneousIC(initial, True)  # cm pressure head, equilibrium
s.setTopBC("noFlux")
s.setBotBC("fleeFlow")
s.setVGParameters([weatherInit['vg']])
s.initializeProblem()
# Sets the critical pressure to limit flow for boundary conditions constantFlow, constantFlowCyl, and atmospheric 
#wilting_point = -15000000  # cm
#s.setCriticalPressure(wilting_point)

sx = s.getSolutionHead()  # inital condition, solverbase.py

picker = lambda x, y, z: s.pick([x, y, z])    
pl.setSoilGrid(picker)  # maps segment
pl2.setSoilGrid(picker)  # maps segment


""" Parameters phloem and photosynthesis """
r = PhloemFluxPython(pl,psiXylInit = min(sx),ciInit = weatherInit["cs"]*0.5) #XylemFluxPython(pl)#
r2 = PhloemFluxPython(pl2,psiXylInit = min(sx),ciInit = weatherInit["cs"]*0.5) #XylemFluxPython(pl)#

setKrKx_phloem()
r.g0 = 8e-3
r.VcmaxrefChl1 =1.28#/2
r.VcmaxrefChl2 = 8.33#/2
r.a1 = 0.5#0#0.6/0.4#0.7/0.3#0.6/0.4 #ci/(cs - ci) for ci = 0.6*cs
r.a3 = 1.5
r.alpha = 0.4#0.2#/2
r.theta = 0.6#0.9#/2
r.k_meso = 1e-3#1e-4
r.setKrm2([[2e-5]])
r.setKrm1([[10e-2]])#([[2.5e-3]])
r.setRhoSucrose([[0.51],[0.65],[0.56]])#0.51
r.setRmax_st([[14.4,9.0,6.0,14.4],[5.,5.],[15.]])#*6 for roots, *1 for stem, *24/14*1.5 for leaves

r.KMrm = 0.1#VERY IMPORTANT TO KEEP IT HIGH
#r.exud_k = np.array([2.4e-4])#*10#*(1e-1)
#r.k_gr = 1#0
r.sameVolume_meso_st = False
r.sameVolume_meso_seg = True
r.withInitVal =True
r.initValST = 0.#0.6#0.0
r.initValMeso = 0.#0.9#0.0
r.beta_loading = 0.6
r.Vmaxloading = 0.05 #mmol/d, needed mean loading rate:  0.3788921068507634
r.Mloading = 0.2
r.Gr_Y = 0.8
r.CSTimin = 0.4
r.surfMeso=0.0025

r.cs = weatherInit["cs"]

#r.r_forPhloem(24/14*1.5, 4)
#r.r_forPhloem(24/14, 3)
#r.r_forPhloem(6, 2) #because roots will have high C_ST less often
r.expression = 6
r.update_viscosity = True
r.solver = 1
r.atol = 1e-12
r.rtol = 1e-8
#r.doNewtonRaphson = False;r.doOldEq = False
SPAD= 41.0
chl_ = (0.114 *(SPAD**2)+ 7.39 *SPAD+ 10.6)/10
r.Chl = np.array( [chl_]) 
r.Csoil = 1e-4


""" for post processing """
structSumInit = 0
orgs_all = r.plant.getOrgans(-1, True)

for org in orgs_all:
    if org.organType() < 2:
        raise Exception("ot < 3")
    structSumInit += org.orgVolume(-1,False) * r.rhoSucrose_f(int(org.getParameter("subType")),
                                                                org.organType())

AnSum = 0
Nt = len(r.plant.nodes) 
Q_Rmbu      = np.array([0.])

Q_Grbu      = np.array([0.])
Q_Gr4bu      = np.array([0.])
Q_Gr3bu      = np.array([0.])
Q_Gr2bu      = np.array([0.])
Q_Grbuth      = np.array([0.])
Q_Gr4buth      = np.array([0.])
Q_Gr3buth      = np.array([0.])
Q_Gr2buth      = np.array([0.])

Q_Exudbu    = np.array([0.])
Q_Rmmaxbu   = np.array([0.])
Q_Grmaxbu   = np.array([0.])
Q_Exudmaxbu = np.array([0.])
Q_STbu      = np.array([0.])
Q_Parbu      = np.array([0.])
Q_mesobu    = np.array([0.])
volSegbu =  np.array([0.])
NOrg = r.plant.getNumberOfOrgans()
delta_ls_bu = np.full(NOrg, 0.)
delta_ls_max = 0


Ntbu = 1
Q_in  = 0
Q_out = 0

#simDuration = 0.

orgs_all = r.plant.getOrgans(-1, True)
volOrgbu = np.array([org.orgVolume(-1,False) for org in orgs_all]) 
volOrg = np.array([org.orgVolume(-1,False) for org in orgs_all]) 
ot_orgs = np.array([org.organType() for org in orgs_all])
st_orgs = np.array([org.getParameter("subType") for org in orgs_all])

volOrgini = np.array([org.orgVolume(-1,False) for org in orgs_all])

volOrgi_th = 0.
lenOrgbu = np.array([org.getLength(False) for org in orgs_all]) 
lenOrg = np.array([org.getLength(False) for org in orgs_all]) 
lenOrgi_th = 0.
Orgidsbu = np.array([org.getId() for org in orgs_all])
Orgids = np.array([org.getId() for org in orgs_all]) #true:realized
#raise Exception
ö=0

orgs_all2 = r2.plant.getOrgans(-1, True)
volOrgini2 = sum(np.array([org.orgVolume(-1 ,True) for org in orgs_all2]) )#true:realized
volOrgi_th2 = 0.


volOrgini2_type =  np.array([sum([org.orgVolume(-1,False) for org in r2.plant.getOrgans(2, True)]),
                             sum([org.orgVolume(-1,False) for org in r2.plant.getOrgans(3, True)]), 
                             sum([org.orgVolume(-1,False) for org in r2.plant.getOrgans(4, True)])]) 
                             
sucOrgini2_type =  np.array([sum([org.orgVolume(-1,False) for org in r2.plant.getOrgans(2, True)])*r.rhoSucrose_f(0,2),
                             sum([org.orgVolume(-1,False) for org in r2.plant.getOrgans(3, True)])*r.rhoSucrose_f(0,3), 
                             sum([org.orgVolume(-1,False) for org in r2.plant.getOrgans(4, True)])*r.rhoSucrose_f(0,4)]) 


                             
volOrgini_type =  np.array([sum([org.orgVolume(-1,False) for org in r.plant.getOrgans(2, True)]),
                            sum([org.orgVolume(-1,False) for org in r.plant.getOrgans(3, True)]), 
                            sum([org.orgVolume(-1,False) for org in r.plant.getOrgans(4, True)])]) 
  
volOrgini2 = np.array([org.orgVolume(-1,False) for org in r.plant.getOrgans(2, True)])
volOrgini3 = np.array([org.orgVolume(-1,False) for org in r.plant.getOrgans(3, True)])
volOrgini4 = np.array([org.orgVolume(-1,False) for org in r.plant.getOrgans(4, True)])
sucOrgini_type =  np.array([sum([org.orgVolume(-1,False) for org in r.plant.getOrgans(2, True)])*r.rhoSucrose_f(0,2),
                            sum([org.orgVolume(-1,False) for org in r.plant.getOrgans(3, True)])*r.rhoSucrose_f(0,3), 
                            sum([org.orgVolume(-1,False) for org in r.plant.getOrgans(4, True)])*r.rhoSucrose_f(0,4)]) 
typeOrg_unit =  np.array([org.organType() for org in r.plant.getOrgans(-1, True)])
idOrg_unit =  np.array([org.getId() for org in r.plant.getOrgans(-1, True)])
sucOrg_unit =  np.array([org.orgVolume(-1,False)*r.rhoSucrose_f(int(org.getParameter("subType")),org.organType()) for org in r.plant.getOrgans(-1, True)])
sucOrgini_unit = sucOrg_unit
#print(volOrgini2_type, volOrgini_type)
sucOrg_type = sucOrgini_type
sucOrg2_type =sucOrgini2_type
volOrg2_type = volOrgini2_type

volOrg_type = volOrgini_type
orgs_roots = r.plant.getOrgans(2, True)
orgs_ln = np.array([])
orgs_st = np.array([])
for org_ in orgs_roots:
    if org_.getParameter("subType") == 2 :
        orgs_ln= np.append(orgs_ln,len(org_.param().ln) ) 
        orgs_st= np.append(orgs_st,org_.getParameter("subType")  ) 

         

orgs_all = r.plant.getOrgans(-1, True)
volOrgbu = np.array([np.full( (len(org.getNodeIds())-1),org.getId())  for org in orgs_all], dtype=object) 

                                                                                               
                                                                                                             
                                                                                                   
                                                                                            

beginning = datetime.now()
#1h for 1d when dxMin = 0.3

AnSum = 0
Q_ST_init = np.array([])
Q_meso_init  = np.array([])
Q_Gr4bu =Q_Gr3bu=Q_Gr2bu=[0]
deltasucorgbu = np.array([])
AnSum = 0
Nt = len(r.plant.nodes) 
Q_Rmbu      = np.array([0.])

Q_Grbu      = np.array([0.])
Q_Gr4bu      = np.array([0.])
Q_Gr3bu      = np.array([0.])
Q_Gr2bu      = np.array([0.])
Q_Grbuth      = np.array([0.])
Q_Gr4buth      = np.array([0.])
Q_Gr3buth      = np.array([0.])
Q_Gr2buth      = np.array([0.])

Q_Exudbu    = np.array([0.])
Q_Rmmaxbu   = np.array([0.])
Q_Grmaxbu   = np.array([0.])
Q_Exudmaxbu = np.array([0.])
Q_STbu      = np.array([0.])
Q_Parbu      = np.array([0.])
Q_mesobu    = np.array([0.])


while simDuration < simMax: 
    
    print('simDuration:',simDuration )
    
    ot_4phloem = r.plant.organTypes # np.insert(,0,2)
    ot_4phloem.insert(0,2)#first node 
    ot_4phloem = np.array(ot_4phloem)
    
    weatherX = weather(simDuration)

    r.Qlight = weatherX["Qlight"] #; TairC = weatherX["TairC"] ; text = "night"
        
        
    setKrKx_xylem(weatherX["TairC"], weatherX["RH"])
    
    r.solve_photosynthesis(sim_time_ = simDuration, sxx_=sx, cells_ = True,RH_ = weatherX["RH"],
        verbose_ = False, doLog_ = False,TairC_= weatherX["TairC"] )
        
    #trans = np.array(r.outputFlux)
    """ dumux """    
    fluxesSoil = r.soilFluxes(simDuration, r.psiXyl, sx, approx=False)
    s.setSource(fluxesSoil.copy())  # richards.py 
    s.solve(dt)
    sx = s.getSolutionHead()  # richards.py    
    min_sx, min_rx, max_sx, max_rx = np.min(sx), np.min(r.psiXyl), np.max(sx), np.max(r.psiXyl)
    n = round((simDuration- simInit)/(simMax-simInit) * 100.)
    print("[" + ''.join(["*"]) * n + ''.join([" "]) * (100 - n) + "], [{:g}, {:g}] cm soil [{:g}, {:g}] cm root at {:g} days {:g}"
            .format(min_sx, max_sx, min_rx, max_rx, s.simTime, r.psiXyl[0]))
          
     
    AnSum += np.sum(r.Ag4Phloem)*dt
    startphloem= simDuration
    endphloem = startphloem + dt
    stepphloem = 1
    filename = "results/pmincpb_" + str(simDuration) + "_15pm.txt" 
    
    errLeuning = sum(r.outputFlux)
    fluxes = np.array(r.outputFlux)
    # save_stdout = sys.stdout
    # sys.stdout = open('trash', 'w')
    verbose_phloem = True
    #deltaVol = r.calcDeltaVolOrgNode(dt,1)
    #print(deltaVol)
    #print("start pm")
    #raise Exception
    filename = "results"+ directoryN +"inPM_"+str(ö)+".txt"
    print("startpm")
    r.startPM(startphloem, endphloem, stepphloem, ( weatherX["TairC"]  +273.15) , verbose_phloem, filename)
    if r.withInitVal and (len(Q_ST_init) ==0) :
        Q_ST_init = np.array(r.Q_init[0:Nt])
        Q_meso_init = np.array(r.Q_init[Nt:(Nt*2)])
        
    #print("start pm_done")
    # sys.stdout = save_stdout
    
    Q_ST    = np.array(r.Q_out[0:Nt])
    Q_meso  = np.array(r.Q_out[Nt:(Nt*2)])
    Q_Rm    = np.array(r.Q_out[(Nt*2):(Nt*3)])
    Q_Exud  = np.array(r.Q_out[(Nt*3):(Nt*4)])
    Q_Gr    = np.array(r.Q_out[(Nt*4):(Nt*5)])
    Q_Gr4       = Q_Gr[np.where(ot_4phloem==4)[0]]#Q_Gr4     - Q_Gr4bu
    Q_Gr3       = Q_Gr[np.where(ot_4phloem==3)[0]]#Q_Gr3     - Q_Gr3bu
    Q_Gr2       = Q_Gr[np.where(ot_4phloem==2)[0]]#Q_Gr2     - Q_Gr2bu
    
    
    C_ST    = np.array(r.C_ST)
    Q_Par   = np.array(r.Q_out[(Nt*8):(Nt*9)])
    Fl      = np.array(r.Fl)
    volST   = np.array(r.vol_ST)
    volMeso   = np.array(r.vol_Meso)
    C_Par   = Q_Par/volST
    C_meso  = Q_meso/volMeso
    Q_in   += sum(np.array(r.AgPhl)*dt)
    Q_out   = Q_Rm + Q_Exud + Q_Gr
    error   = sum(Q_ST +Q_Par+ Q_meso + Q_out )- Q_in - sum(Q_ST_init)  - sum(Q_meso_init)
    # Q_ST_dot    = np.array(r.Q_out_dot[0:Nt])
    # Q_meso_dot  = np.array(r.Q_out_dot[Nt:(Nt*2)])
    # Q_Rm_dot    = np.array(r.Q_out_dot[(Nt*2):(Nt*3)])
    # Q_Exud_dot  = np.array(r.Q_out_dot[(Nt*3):(Nt*4)])
    # Q_Gr_dot    = np.array(r.Q_out_dot[(Nt*4):(Nt*5)])
    
    #delta_ls_max += sum(np.array(r2.rmaxSeg(dt, r.k_gr)) * dt)
    delta_ls_max_i = np.array(r.delta_ls_org_imax)
    delta_ls_max = np.array(r.delta_ls_org_max)
    delta_ls_i = np.array(r.delta_ls_org_i)
    delta_ls = np.array(r.delta_ls_org)
    
    Q_Rmmax       = np.array(r.Q_out[(Nt*5):(Nt*6)])
    Q_Grmax       = np.array(r.Q_out[(Nt*6):(Nt*7)])
    Q_Exudmax     = np.array(r.Q_out[(Nt*7):(Nt*8)])
    
    Q_ST_i        = Q_ST      - Q_STbu
    Q_Par_i       = Q_out     - Q_Parbu
    Q_Rm_i        = Q_Rm      - Q_Rmbu
    Q_Gr_i        = Q_Gr      - Q_Grbu
    Q_Gr4_i       = Q_Gr_i[np.where(ot_4phloem==4)[0]]#Q_Gr4     - Q_Gr4bu
    Q_Gr3_i       = Q_Gr_i[np.where(ot_4phloem==3)[0]]#Q_Gr3     - Q_Gr3bu
    Q_Gr2_i       = Q_Gr_i[np.where(ot_4phloem==2)[0]]#Q_Gr2     - Q_Gr2bu
    
    #Q_Gr_ith        = Q_Grth      - Q_Grbuth
    #Q_Gr4_ith       = Q_Gr4th     - Q_Gr4buth
    #Q_Gr3_ith       = Q_Gr3th     - Q_Gr3buth
    #Q_Gr2_ith       = Q_Gr2th     - Q_Gr2buth
    
    Q_Exud_i      = Q_Exud    - Q_Exudbu
    Q_meso_i      = Q_meso    - Q_mesobu
    
    Q_Rmmax_i     = Q_Rmmax   - Q_Rmmaxbu
    Q_Grmax_i     = Q_Grmax   - Q_Grmaxbu
    Q_Exudmax_i   = Q_Exudmax - Q_Exudmaxbu
    
    Q_out_i       = Q_Rm_i    + Q_Exud_i      + Q_Gr_i
    Q_outmax_i    = Q_Rmmax_i + Q_Exudmax_i   + Q_Grmax_i
    
    
    orgs = r.plant.getOrgans(-1, True)
    id_orgs = np.array([org.getId() for org in orgs])
    orgs_all = r.plant.getOrgans(-1, True)
    ot_orgs_all = np.array([org.organType() for org in orgs_all])
    volOrgi_th = 0#volOrg2 - volOrgini2
    
    volOrg2 = np.array([org.orgVolume(-1,False) for org in r.plant.getOrgans(2, True)])
    volOrg3 = np.array([org.orgVolume(-1,False) for org in r.plant.getOrgans(3, True)])
    volOrg4 = np.array([org.orgVolume(-1,False) for org in r.plant.getOrgans(4, True)])
    
    volOrg_typei = volOrg_type - volOrgini_type
    volOrg2_typei = volOrg2_type - volOrgini2_type
    
    JW_ST = np.array(r.JW_ST)
    length_ST = np.array(r.segLength())
    #0.0001037
    Lp = 0.005#0.004320 #cm d-1 hPa, assumed resistance between xylem and phloem
    Rgaz =  83.14 #hPa cm3 K-1 mmol-1
    a_STs = np.array(r.a_ST)#np.array([a_ST[ot][st] for ot, st in ])
    #RhatFhat =   (weatherX["TairC"] + 273.15) * C_ST[1:] * Rgaz * (2/a_STs) * length_ST * Lp /np.abs(JW_ST) /   ((a_STs**2)*np.pi)  
    #first part of RhatFhat: 9compute rest afterwards
    RhatFhat =   (weatherX["TairC"] + 273.15) * C_ST[1:] * Rgaz * length_ST  /np.abs(JW_ST) #* Lp *  perimeter  
    C_ST_ = C_ST[1:]
    ids = np.where(RhatFhat ==  min(RhatFhat))
    if (min(RhatFhat) < 1) :
        #print()
        C_ST_ = C_ST[1:]
        ids = np.where(RhatFhat ==  min(RhatFhat))
        print(min(RhatFhat))
        print(C_ST_[ids] , Rgaz  )
        print(a_STs[ids]  )
        print(length_ST[ids]   )
        print(JW_ST[ids]) 
        print( RhatFhat[ids],(weatherX["TairC"] + 273.15)  )
        print("issue RhatFhat")
        #raise Exception
    
    if verbose :
        print("\n\n\n\t\tat ", int(np.floor(simDuration)),"d", int((simDuration%1)*24),"h",  round(r.Qlight *1e6),"mumol m-2 s-1")
        print("Error in Suc_balance:\n\tabs (mmol) {:5.2e}\trel (-) {:5.2e}".format(error, div0f(error,Q_in, 1.)))
        #print("Error in growth:\n\tabs (mmol) {:5.2e}\trel (-) {:5.2e}".format(errorGri, relErrorGri))
        print("Error in photos:\n\tabs (cm3/day) {:5.2e}".format(errLeuning))
        print("water fluxes (cm3/day):\n\ttrans {:5.2e}\tminExud {:5.2e}\tmaxExud {:5.2e}".format(sum(fluxesSoil.values()), min(fluxesSoil.values()), max(fluxesSoil.values())))
        print("C_ST (mmol ml-1):\n\tmean {:.2e}\tmin  {:5.2e} at {:d} segs \tmax  {:5.2e}".format(np.mean(C_ST), min(C_ST), len(np.where(C_ST == min(C_ST) )[0]), max(C_ST)))        
        print("C_me (mmol ml-1):\n\tmean {:.2e}\tmin  {:5.2e}\tmax  {:5.2e}".format(np.mean(C_meso), min(C_meso), max(C_meso)))        
        print('Q_X (mmol Suc): \n\tST   {:.2e}\tmeso {:5.2e}\tin   {:5.2e}'.format(sum(Q_ST), sum(Q_meso), Q_in))
        print('\tRm   {:.2e}\tGr   {:.2e}\tExud {:5.2e}'.format(sum(Q_Rm), sum(Q_Gr), sum(Q_Exud)))
        print("aggregated sink satisfaction at last time step (%) :\n\ttot  {:5.1f}\n\tRm   {:5.1f}\tGr   {:5.1f}\tExud {:5.1f}".format(
            sum(Q_out_i)/sum(Q_outmax_i)*100,sum(Q_Rm_i)/sum(Q_Rmmax_i)*100, 
             div0f(sum(Q_Gr_i), sum(Q_Grmax_i), 1.)*100,div0f(sum(Q_Exud_i),sum(Q_Exudmax_i), 1.)*100))
        print("aggregated sink repartition at last time step (%) :\n\tRm   {:5.1f}\tGr   {:5.1f}\tExud {:5.1f}".format(sum(Q_Rm_i)/sum(Q_out_i)*100, 
             sum(Q_Gr_i)/sum(Q_out_i)*100,sum(Q_Exud_i)/sum(Q_out_i)*100))
        print("aggregated sink repartition (%) :\n\tRm   {:5.1f}\tGr   {:5.1f}\tExud {:5.1f}".format(sum(Q_Rm)/sum(Q_out)*100, 
             sum(Q_Gr)/sum(Q_out)*100,sum(Q_Exud)/sum(Q_out)*100))
        print("aggregated sink repartition for max (%) :\n\tRm   {:5.1f}\tGr   {:5.1f}\tExud {:5.1f}".format(sum(Q_Rmmax_i)/sum(Q_outmax_i)*100, 
             sum(Q_Grmax_i)/sum(Q_outmax_i)*100,sum(Q_Exudmax_i)/sum(Q_outmax_i)*100))
        print("abs val for max :\n\tRm   {:5.5f}\tGr   {:5.5f}\tExud {:5.5f}".format(sum(Q_Rmmax_i), 
             sum(Q_Grmax_i),sum(Q_Exudmax_i)))
        print("Q_Par {:5.2e}, C_Par {:5.2e}".format(sum(Q_Par), np.mean(C_Par)))
        print("growth (cm)\ttot {:5.2e}\ti {:5.2e}".format(sum(delta_ls), sum(delta_ls_i)))      
        print("max growth (cm)\ttot {:5.2e}\ti {:5.2e}".format(sum(delta_ls_max), sum(delta_ls_max_i)))     
        print("amount Suc (cm)\tAn {:5.2e}\tGr {:5.2e}\tRGr {:5.2e}\tRm {:5.2e}\tExud {:5.2e}".format(AnSum, sum(Q_Gr)*r.Gr_Y,sum(Q_Gr)*(1-r.Gr_Y), sum(Q_Rm), sum(Q_Exud))) 
        #print("growth (cm3)\n\tobs {:5.2e}\ttotth {:5.2e}\ttotobs {:5.2e}".format(sum(volOrgi_obs), volOrgi_th,volOrg_obstot ))  
        print("growth (cm3) per type\n\ttotobs", volOrg_typei , volOrg2_typei)       
        print("sucOrg obs (mmol)\t th (mmol)\t", sucOrg_type - sucOrgini_type, sucOrg2_type - sucOrgini2_type)       
        print("Grobs (mmol) root\tstem\tleaf\t", sum(Q_Gr2bu)*r.Gr_Y,sum(Q_Gr3bu)*r.Gr_Y, sum(Q_Gr4bu)*r.Gr_Y)# , gr4ith) 
        #print("RhatFhat ", min(RhatFhat),C_ST_[ids],a_STs[ids], length_ST[ids], JW_ST[ids]  )
        
    if min(C_ST) < 0.0:
        print("min(C_ST) < 0.015", min(C_ST),np.mean(C_ST),max(C_ST))
        raise Exception
    write_file_array("RhatFhat", RhatFhat)
    write_file_array("a_STs", a_STs)
    write_file_array("JW_ST", JW_ST)#cm3/d
    write_file_array("length_ST", length_ST)
    #NB: delta_ls_max inexact as does not account for growth of organxs which are not here at the beginning
    #print([org.orgVolume(-1,False) for org in r.plant.getOrgans(3, True)],[org.orgVolume(-1,False) for org in r2.plant.getOrgans(3, True)])
    #print([org.getParameter("lengthTh") for org in r.plant.getOrgans(3, True)],[org.getParameter("lengthTh") for org in r2.plant.getOrgans(3, True)])
    # print([org.getId() for org in r.plant.getOrgans(3, True)],[org.orgVolume(-1,False) for org in r2.plant.getOrgans(3, True)])
    # print(len(r.plant.getOrgans(2, True)),len(r.plant.getOrgans(3, True)),len(r.plant.getOrgans(4, True)))
    #dist2tip = np.array(r.plant.dist2tips)
    ana = pb.SegmentAnalyser(r.plant.mappedSegments())
    
    cutoff = 1e-15 #is get value too small, makes paraview crash
    C_ST_p = C_ST
    C_ST_p[abs(C_ST_p) < cutoff] = 0
    fluxes_p = fluxes
    fluxes_p[abs(fluxes_p) < cutoff] = 0
    Q_Exud_i_p = Q_Exud_i
    Q_Exud_i_p[abs(Q_Exud_i_p) < cutoff] = 0
    Q_Rm_i_p = Q_Rm_i
    Q_Rm_i_p[abs(Q_Rm_i_p) < cutoff] = 0
    Q_Gr_i_p = Q_Gr_i
    Q_Gr_i_p[abs(Q_Gr_i_p) < cutoff] = 0
    
    Q_Exudmax_i_p = Q_Exudmax_i
    Q_Exudmax_i_p[abs(Q_Exudmax_i_p) < cutoff] = 0
    Q_Rmmax_i_p = Q_Rmmax_i
    Q_Rmmax_i_p[abs(Q_Rmmax_i_p) < cutoff] = 0
    Q_Grmax_i_p = Q_Grmax_i
    Q_Grmax_i_p[abs(Q_Grmax_i_p) < cutoff] = 0
    
    
    C_Exud_i_p = Q_Exud_i/volST
    C_Exud_i_p[abs(C_Exud_i_p ) < cutoff] = 0
    C_Rm_i_p = Q_Rm_i/volST
    C_Rm_i_p[abs(C_Rm_i_p) < cutoff] = 0
    C_Gr_i_p = Q_Gr_i/volST
    C_Gr_i_p[abs(C_Gr_i_p) < cutoff] = 0
    
    C_Exudmax_i_p = Q_Exudmax_i/volST
    C_Exudmax_i_p[abs(C_Exudmax_i_p) < cutoff] = 0
    C_Rmmax_i_p = Q_Rmmax_i/volST
    C_Rmmax_i_p[abs(C_Rmmax_i_p) < cutoff] = 0
    C_Grmax_i_p = Q_Grmax_i/volST
    C_Grmax_i_p[abs(C_Grmax_i_p) < cutoff] = 0
    
    psiXyl_p = np.array(r.psiXyl)
    psiXyl_p[abs(psiXyl_p) < cutoff] = 0
    ana.addData("CST", C_ST_p)
    #do as konz or div per vol or surf?
    #ana.addData("Q_Exud", Q_Exud)  # cut off for vizualisation
    ana.addData("fluxes", fluxes_p)
    ana.addData("Fpsi", np.array(r.Fpsi))
    
    ana.addData("QExud", Q_Exud_i_p)  # cut off for vizualisation
    ana.addData("QRm", Q_Rm_i_p)  # cut off for vizualisation
    ana.addData("QGr", Q_Gr_i_p)  # cut off for vizualisation
    ana.addData("QExudmax", Q_Exudmax_i_p)  # cut off for vizualisation
    ana.addData("QRmmax", Q_Rmmax_i_p)  # cut off for vizualisation
    ana.addData("QGrmax", Q_Grmax_i_p)  # cut off for vizualisation
    
    ana.addData("CExud", C_Exud_i_p)  # cut off for vizualisation
    ana.addData("CRm", C_Rm_i_p)  # cut off for vizualisation
    ana.addData("CGr", C_Gr_i_p)  # cut off for vizualisation
    ana.addData("CExudmax", C_Exudmax_i_p)  # cut off for vizualisation
    ana.addData("CRmmax", C_Rmmax_i_p)  # cut off for vizualisation
    ana.addData("CGrmax", C_Grmax_i_p)  # cut off for vizualisation
    
    ana.addData("psi_Xyl",psiXyl_p)
    ana.write("results"+directoryN+"15pm_"+ str(ö) +".vtp", ["CST", "fluxes","psi_Xyl",
                        "QExud", "QGr", "QRm",
                        "CExud", "CGr", "CRm",
                        "QExudmax", "QGrmax", "QRmmax",
                        "CExudmax", "CGrmax", "CRmmax",
                        "organType", "subType", "Fpsi"]) 
    
    rl_ = ana.distribution("length", 0., -depth, layers, True)                   
    rl_ = np.array(rl_)/ soilvolume  # convert to density
    fluxes_p = np.insert(fluxes_p,0,0)# "[sucrose]",
    if(simDuration > (simMax -1)):
        vp.plot_plant(r.plant,p_name = [ "fluxes","Exud","Gr","Rm","xylem pressure (cm)"],
                            vals =[ fluxes_p, Q_Exud_i_p, Q_Gr_i_p, Q_Rm_i_p, psiXyl_p], 
                            filename = "results"+ directoryN +"plotplant_psi_"+ str(ö), range_ = [300,1450])
        vp.plot_plant(r.plant,p_name = [ "fluxes","Exud","Gr","Rm","sucrose concentration (mmol/cm3)"],
                            vals =[ fluxes_p, Q_Exud_i_p, Q_Gr_i_p, Q_Rm_i_p, C_ST_p], 
                            filename = "results"+ directoryN +"plotplant_suc_"+ str(ö), range_ = [0,3])   
    ö +=1
    r_ST_ref = np.array(r.r_ST_ref)
    r_ST = np.array(r.r_ST)
    
    #with open("results"+ directoryN +"CWGr_max_15pm.txt",  'a') as data: 
      #    data.write(str(r.deltaSucOrgNode))
          
    ots = np.concatenate((np.array([0]), r.get_organ_types()))#per node
    #ots_org = np.concatenate((np.array([0]), r.get_organ_types()))#per org
    #write_file_array("id_org", ot_orgs_all)
    #write_file_array("ots_org", ot_orgs_all)
    #sucOrg_unit =  np.array([org.orgVolume(-1,False)*r.rhoSucrose[org.organType()] for org in r.plant.getOrgans(-1, True)])
    #typeOrg_unit =  np.array([org.organType() for org in r.plant.getOrgans(-1, True)])
    
    write_file_array("RLD", rl_)
    write_file_array("Fpsi", r.Fpsi)
    write_file_array("deltasucorgbu_type", typeOrg_unit)
    write_file_array("deltasucorgbu_phloem", deltasucorgbu)
    write_file_array("deltasucorgbu_plant",  (sucOrg_unit - sucOrgini_unit)[idOrg_unit.argsort()])
    write_file_array("deltasucorgbu_plantid",  idOrg_unit[idOrg_unit.argsort()])
    
    #write_file_array("volOrg2", volOrg2-volOrgini2)
    #write_file_array("volOrg3", volOrg3-volOrgini3)
    #write_file_array("volOrg4", volOrg4-volOrgini4)
    
    write_file_array("leafBladeSurface", np.array(r.plant.leafBladeSurface))
    write_file_array("fw", r.fw)
    write_file_array("gco2", r.gco2)
    
    #write_file_array("length_blade", length_blade)
    write_file_array("ots", ots)
    write_file_array("soilWatPot", sx)
    write_file_array("fluxes", fluxes)#cm3 day-1
    write_file_array("volST", volST)
    write_file_array("volOrg",  volOrg) #with epsilon
    write_file_array("Fl", Fl)
    write_file_array("AgPhl", np.array(r.AgPhl))
    write_file_array("Q_ST_dot", Q_ST_i/dt)
    write_file_array("Q_meso_dot", Q_meso_i/dt)
    write_file_array("Q_Rm_dot", Q_Rm_i/dt)
    write_file_array("Q_Exud_dot", Q_Exud_i/dt)
    write_file_array("Q_Gr_dot", Q_Gr_i/dt)
    write_file_array("Q_Rmmax_dot", Q_Rmmax_i/dt)
    write_file_array("Q_Exudmax_dot", Q_Exudmax_i/dt)
    write_file_array("Q_Grmax_dot", Q_Grmax_i/dt)
    
    write_file_array("Q_Par", Q_Par)
    write_file_array("C_Par", C_Par)
    write_file_array("Q_ST", Q_ST)
    write_file_array("C_ST", C_ST)
    write_file_array("C_meso", C_meso)
    write_file_array("Q_meso", Q_meso)
    write_file_array("Q_Rm", Q_Rm)
    write_file_array("Q_Exud", Q_Exud)
    write_file_array("Q_Gr", Q_Gr)
    write_file_array("psiXyl", r.psiXyl)
    write_file_array("trans", r.Ev)
    write_file_array("transrate",r.Jw)
    leavesSegs = np.where(ots[1:] ==4)
    fluxes_leaves = fluxes[leavesSegs]
    if (min(r.Ev) < 0) or (min(r.Jw) < 0) or (min(fluxes_leaves)<0):
        print("leaf looses water", min(r.Ev),min(r.Jw), min(fluxes_leaves))
        raise Exception
    write_file_array("Q_Grmax", Q_Grmax)
    write_file_array("Q_Rmmax", Q_Rmmax)
    write_file_array("Q_Exudmax", r.Q_Exudmax)
    
    
    write_file_array("ratio_Gr ", div0(Q_Gr_i,Q_Grmax_i, np.nan))
    write_file_array("ratio_Rm ", Q_Rm_i/Q_Rmmax_i)
    write_file_array("ratio_Exud ", div0(Q_Exud_i,Q_Exudmax_i, np.nan))
    write_file_array("satisfaction ", Q_out_i/Q_outmax_i)
    
    write_file_array("r_ST_ref", r_ST_ref)
    write_file_array("r_ST", r_ST)
    write_file_array("mu", r_ST/r_ST_ref)#hPa d
    
    write_file_array("id_orgs", id_orgs)
    write_file_array("ot_orgs", ot_orgs)
    write_file_array("ot_orgs_all", ot_orgs)#with arg to small to be represented
    write_file_array("st_orgs", st_orgs)
    write_file_array("delta_ls", delta_ls)
    write_file_array("Q_Ag", r.AgPhl)
    write_file_array("delta_ls_max", delta_ls_max)
    write_file_array("delta_ls_i", delta_ls_i)
    write_file_array("delta_ls_max_i", delta_ls_max_i)
    write_file_float("time", simDuration)
    write_file_float("computeTime", datetime.now() - beginning)
    organTypes = r.get_organ_types()
    subTypes =r.get_subtypes()
    nodes_organtype = np.concatenate(([1], organTypes))#seg_Indx = node_y_Indx - 1. just need to add organ type for seed (set to type 'root')     #np.zeros(len(nods))#-1)
    nodes_subtype = np.concatenate(([1], subTypes))
    write_file_array("organTypes ", organTypes)
    write_file_array("subTypes ", subTypes)
    
    Ntbu = Nt
    #Ntbu2 = Nt2
    orgs_all = r.plant.getOrgans(-1, True)
    NOrgbu = len(orgs_all)
    orgradibu= np.array([org.getParameter("a") for org in orgs_all])
    volOrgbu = np.array([org.orgVolume(-1,False) for org in orgs_all])
    lenOrgbu = np.array([org.getLength(False) for org in orgs_all])
    Orgidsbu = np.array([org.getId() for org in orgs_all]) 
    
    verbose_simulate = False
    r.plant.simulate(dt, verbose_simulate)#, "outputpm15.txt") #time span in days /(60*60*24)
    r2.plant.simulate(dt,  verbose_simulate)#, "outputpm15.txt")
    
    volOrg2_type =  np.array([sum([org.orgVolume(-1,False) for org in r2.plant.getOrgans(2, True)]),
                            sum([org.orgVolume(-1,False) for org in r2.plant.getOrgans(3, True)]), 
                            sum([org.orgVolume(-1,False) for org in r2.plant.getOrgans(4, True)])]) 
    
    
    volOrg_type =  np.array([sum([org.orgVolume(-1,False) for org in r.plant.getOrgans(2, True)]),
                            sum([org.orgVolume(-1,False) for org in r.plant.getOrgans(3, True)]), 
                            sum([org.orgVolume(-1,False) for org in r.plant.getOrgans(4, True)])])
    
    
    # lenOrg_type =  np.array([sum([org.orgVolume(-1,False) for org in r.plant.getOrgans(2, True)]),
                            # sum([org.orgVolume(-1,False) for org in r.plant.getOrgans(3, True)]), 
                            # sum([org.orgVolume(-1,False) for org in r.plant.getOrgans(4, True)])]) 
    
    sucOrg2_type =  np.array([sum([org.orgVolume(-1,False) for org in r2.plant.getOrgans(2, True)])*r.rhoSucrose_f(0,2),
                            sum([org.orgVolume(-1,False) for org in r2.plant.getOrgans(3, True)])*r.rhoSucrose_f(0,3), 
                            sum([org.orgVolume(-1,False) for org in r2.plant.getOrgans(4, True)])*r.rhoSucrose_f(0,4)]) 
    sucOrg_type =  np.array([sum([org.orgVolume(-1,False) for org in r.plant.getOrgans(2, True)])*r.rhoSucrose_f(0,2),
                            sum([org.orgVolume(-1,False) for org in r.plant.getOrgans(3, True)])*r.rhoSucrose_f(0,3), 
                            sum([org.orgVolume(-1,False) for org in r.plant.getOrgans(4, True)])*r.rhoSucrose_f(0,4)]) 
                            
    sucOrg_unit =  np.array([org.orgVolume(-1,False)*r.rhoSucrose_f(int(org.getParameter("subType")),org.organType()) for org in r.plant.getOrgans(-1, True)])
    typeOrg_unit =  np.array([org.organType() for org in r.plant.getOrgans(-1, True)])
    idOrg_unit =  np.array([org.getId() for org in r.plant.getOrgans(-1, True)])
    
    deltasucorgbu = np.array(r.delta_suc_org)   
    
    write_file_array("volOrgth", volOrg2_type)
    write_file_array("volOrgobs", volOrg_type)
    write_file_array("sucOrgth", sucOrg2_type)
    write_file_array("sucOrgobs", sucOrg_type)
    
    orgs_all2 = r2.plant.getOrgans(-1, True)
    volOrgi_th2 =  sum(np.array([org.orgVolume(-1,False) for org in orgs_all2]) )#true:realized
    
    orgs = r.plant.getOrgans(-1, True)
    orgs_all = r.plant.getOrgans(-1, True)
    
    Nt = len(r.plant.nodes)
    NOrg = len(orgs_all)#r.plant.getNumberOfOrgans()#att! not same number
    sucOrgini_unit = np.concatenate((sucOrgini_unit, np.full(NOrg - NOrgbu, 0.)))
    
    orgradi= np.array([org.getParameter("a") for org in orgs_all])
    orgradibu = np.concatenate((orgradibu, np.full(NOrg - NOrgbu, 0.)))
    #Orgidsbu = Orgids
    Orgids = np.array([org.getId() for org in orgs_all]) 
    
    Orgidsbu = np.concatenate((Orgidsbu, np.full(NOrg - NOrgbu, 0.)))
    
    volOrg = np.array([org.orgVolume(-1,False) for org in orgs_all])
    ot_orgs = np.array([org.organType() for org in orgs_all])
    st_orgs = np.array([org.getParameter("subType") for org in orgs_all])
    stroot = st_orgs[np.where(ot_orgs==2)]
    # errorst = sum(np.where(stroot > 3))
    # if errorst > 0:
        # raise Exception
    
    volOrgbu     = np.concatenate((volOrgbu, np.full(NOrg - NOrgbu, 0.)))#len(volOrg) - len(volOrgbu), 0.)))
    #volOrgi_th = sum(Q_Gr_i)*r.Gr_Y/r.rhoSucrose
    
    
    #not sur I get the organs always in same order
    lenOrg = np.array([org.getLength(False) for org in orgs_all]) #true:realized
    lenOrgbu     = np.concatenate((lenOrgbu, np.full(NOrg - NOrgbu, 0.)))#len(volOrg) - len(volOrgbu), 0.)))
    lenOrgi_th = np.concatenate((delta_ls_i, np.full(NOrg - NOrgbu, 0.)))
    
    #delta_ls_max =  np.concatenate((delta_ls_max, np.full(Nt2 - Ntbu2, 0.)))
    Q_Rmbu       =   np.concatenate((Q_Rm, np.full(Nt - Ntbu, 0.)))
    Q_Grbu       =   np.concatenate((Q_Gr, np.full(Nt - Ntbu, 0.))) 
    Q_Exudbu     =   np.concatenate((Q_Exud, np.full(Nt - Ntbu, 0.))) 
    
    Q_Rmmaxbu    =   np.concatenate((Q_Rmmax, np.full(Nt - Ntbu, 0.)))
    Q_Grmaxbu    =   np.concatenate((Q_Grmax, np.full(Nt - Ntbu, 0.))) 
    Q_Exudmaxbu  =   np.concatenate((Q_Exudmax, np.full(Nt - Ntbu, 0.))) 
    
    Q_STbu       =   np.concatenate((Q_ST, np.full(Nt - Ntbu, 0.)))
    Q_mesobu     =   np.concatenate((Q_meso, np.full(Nt - Ntbu, 0.)))
    delta_ls_bu  =   np.concatenate((delta_ls, np.full(NOrg - NOrgbu, 0.)))
    
    Q_Gr4bu= Q_Gr4
    Q_Gr3bu =Q_Gr3
    Q_Gr2bu= Q_Gr2
    ##
    #       CHECKS
    ##
    highNeed1 = sum(Q_Rm_i[np.greater(Q_Rm_i,Q_Rmmax_i)] - Q_Rmmax_i[np.greater(Q_Rm_i,Q_Rmmax_i)])
    highNeed2 = sum(Q_Gr_i[np.greater(Q_Gr_i,Q_Grmax_i)] - Q_Grmax_i[np.greater(Q_Gr_i,Q_Grmax_i)])
    highNeed3 = sum(Q_Exud_i[np.greater(Q_Exud_i,Q_Exudmax_i)] - Q_Exudmax_i[np.greater(Q_Exud_i,Q_Exudmax_i)])
    
    if highNeed1 > abs(error):
        print(np.where([np.greater(Q_Rm_i,Q_Rmmax_i)]))
    #assert div0f(error,Q_in,1.) < 1e-3, "balance error > 0.1%"
    
    if (len(Orgids) != len(Orgidsbu)):
        print(len(Orgids), len(Orgidsbu), NOrg, NOrgbu, Nt, Ntbu)
    orderId  = Orgids - Orgidsbu
    
    lenOrgi_obs = lenOrg - lenOrgbu
    lenOrgi_obs2 = sum(lenOrg) - sum(lenOrgbu)
    errorleni = lenOrgi_obs -  lenOrgi_th 
    errorleni2 = lenOrgi_obs2 -  sum(lenOrgi_th) 
    
    simDuration += dt
    #raise Exception
    # if  abs(errorleni2) > 1e-3:
        # print(lenOrgi_obs)
        # print(lenOrgi_th)
        # print(errorleni)
        # print(lenOrgi_obs2,sum(lenOrgi_th) )
        # print(errorleni2)
        # print( "len error")
        # raise Exception
    
    
    # volOrgi_obs = volOrg - volOrgbu
    # volOrgi_obs2 = sum(volOrg) - sum(volOrgbu)
    # errorGri = volOrgi_obs2 -  volOrgi_th 
    # relErrorGri = div0f(errorGri,volOrgi_th ,1.)
    # assert abs(errorGri) < 1e-3, "vol errorend"
structSum = 0
orgs_all = r.plant.getOrgans(-1, True)

for org in orgs_all: #cm3 * mmol Suc /cm3 = mmol Suc 
    structSum += org.orgVolume(-1,False) * r.rhoSucrose_f(int(org.getParameter("subType")),org.organType())


GrowthSum = structSum - structSumInit
print("simDuration", simDuration, "d")
#print("AnSum", AnSum, "mmol Suc")
#print("GrowthSum", GrowthSum*(1/0.75),"resp ",AnSum*0.95 - GrowthSum*(1/0.8), (AnSum*0.95 - GrowthSum*(1/0.75))/AnSum*100, "exud ", AnSum*0.05, "mmol Suc", "init suc",structSumInit, "end sum", structSum)

end = datetime.now()
print(end - beginning)