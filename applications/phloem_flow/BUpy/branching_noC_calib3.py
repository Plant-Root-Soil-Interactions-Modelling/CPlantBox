""" water movement within the root (static soil) """


import sys; 
#directoryN = "/"+sys.argv[0].split('.')[0]+"/"
directoryN = "/branching_noC/"
sys.path.append("../.."); sys.path.append("../../src/python_modules")
CPBdir = "../.."
#import matplotlib
#matplotlib.use('AGG') #otherwise "import matplotlib.pyplot" hangs

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
        if item.startswith("/results"+directoryN+"prints"):
            os.remove(os.path.join(dir_name2, item))
delete_verbose()


isCluster = (os.environ['HOME'] == '/home/m.giraud')
if False:
    def print(*args, **kwargs):
        """ custom print() function.
            for cluster: can get output even if program stop
            unexpectedly (i.e., without creating the outputfile)
        """
        # Adding new arguments to the print function signature
        # is probably a bad idea.
        # Instead consider testing if custom argument keywords
        # are present in kwargs
        if 'sep' in kwargs:
            sep = kwargs['sep']
        else:
            sep = ' '
        home_dir = os.getcwd()
        dir_name =  "/results"+directoryN
        dir_name2 = home_dir + dir_name
        name2 = dir_name2 + 'prints.txt'
        with open(name2, 'a') as log:
            for arg in args: log.write(str(arg) + sep)
            log.write('\n')

#time lapse pea plants:
#https://www.youtube.com/watch?v=a2Mxpg5lNFM
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
def stairs(t):
    hours = t%1
    coef = int(hours>=0.5)
    return coef

def qair2rh(qair, es_, press):
    e =qair * press / (0.378 * qair + 0.622)
    rh = e / es_
    rh=max(min(rh, 1.),0.)
    return rh


def weather(simDuration):
    vgSoil = [0.059, 0.45, 0.00644, 1.503, 1]
    
    Qmin = 0; Qmax = 960e-6 #458*2.1
    Tmin = 22; Tmax = 25
    specificHumidity = 0.0097
    Pair = 1010.00 #hPa
    thetaInit = 30/100

    coefhours = stairs(simDuration)#sinusoidal(simDuration)
    TairC_ = Tmin + (Tmax - Tmin) * coefhours
    Q_ = Qmin + (Qmax - Qmin) * coefhours
    cs = 350e-6 #co2 paartial pressure at leaf surface (mol mol-1)
    #RH = 0.5 # relative humidity
    es =  6.112 * np.exp((17.67 * TairC_)/(TairC_ + 243.5))
    RH = 0.6#qair2rh(specificHumidity, es, Pair)
    
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
    l_kr = 100 #cm
    #r.setKr([[kr_r0], [kr_s], [kr_l]], kr_length_=l_kr)
    #r.setKx([[kz_r0],[kz_s],[kz_l]])
    r.setKr([[kr_r0,kr_r1,kr_r2,kr_r0],[kr_s,kr_s ,kr_s,kr_s ],[kr_s,kr_l]], kr_length_=l_kr)
    r.setKx([[kz_r0,kz_r1,kz_r2,kz_r0],[kz_s,kz_s,kz_s,kz_s ],[kz_s,kz_l]])
    
    
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
simInit = 10
simDuration = simInit # [day] init simtime
simMax =11
depth = 60
dt = 1/24 #1h
verbose = True

# plant system 
pl = pb.MappedPlant(seednum = 2) #pb.MappedRootSystem() #pb.MappedPlant()
pl2 = pb.MappedPlant(seednum = 2) #pb.MappedRootSystem() #pb.MappedPlant()
path = CPBdir+"/modelparameter/plant/"
name = "morning_glory_UQ"#"wheat_uqr15" #"manyleaves"## "Anagallis_femina_Leitner_2010"  # Zea_mays_1_Leitner_2010

pl.readParameters(path + name + ".xml")
pl2.readParameters(path + name + ".xml")


#raise Exception
sdf = pb.SDF_PlantBox(np.Inf, np.Inf, depth )

pl.setGeometry(sdf) # creates soil space to stop roots from growing out of the soil
pl2.setGeometry(sdf) # creates soil space to stop roots from growing out of the soil


pl.initialize(verbose = False)#, stochastic = False)
pl.simulate(simDuration, False)#, "outputpm15.txt")
pl2.initialize(verbose = False)#, stochastic = False)
pl2.simulate(simDuration, False)#, "outputpm15.txt")

ana = pb.SegmentAnalyser(pl.mappedSegments())
ana.write("results"+directoryN+name + "_"+ str(0) +".vtp") 

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

# s = RichardsWrapper(RichardsSP())
# s.initialize()
# periodic = True
# s.createGrid(min_b, max_b, cell_number, periodic)  # [cm]
# s.setHomogeneousIC(initial, True)  # cm pressure head, equilibrium
# s.setTopBC("noFlux")
# s.setBotBC("fleeFlow")
# s.setVGParameters([weatherInit['vg']])
# s.initializeProblem()
# # Sets the critical pressure to limit flow for boundary conditions constantFlow, constantFlowCyl, and atmospheric 
# #wilting_point = -15000000  # cm
# #s.setCriticalPressure(wilting_point)

# sx = s.getSolutionHead()  # inital condition, solverbase.py
p_mean = weatherInit["p_mean"]#-187
p_bot = p_mean + depth/2
p_top = p_mean - depth/2
sx = np.linspace(p_top, p_bot, depth)
print(sx)
#p_g = -2000 # water potential of the gua
picker = lambda x,y,z : max(int(np.floor(-z)),-1) #abovegroud nodes get index -1

#picker = lambda x, y, z: s.pick([x, y, z])    
pl.setSoilGrid(picker)  # maps segment
pl2.setSoilGrid(picker)  # maps segment


""" Parameters phloem and photosynthesis """
r = PhloemFluxPython(pl,psiXylInit = min(sx),ciInit = weatherInit["cs"]*0.5) #XylemFluxPython(pl)#
r2 = PhloemFluxPython(pl2,psiXylInit = min(sx),ciInit = weatherInit["cs"]*0.5) #XylemFluxPython(pl)#

setKrKx_phloem()
r.k_meso = 1e-3#1e-4
r.setKrm2([[2e-5]])
r.setKrm1([[10e-2]])#([[2.5e-2]])
r.setRhoSucrose([[0.1],[0.1],[0.1]])#0.51
r.setRmax_st([[4.,2.,1.0,4.],[4.,0.00000001],[2.,2.]])#*6 for roots, *1 for stem, *24/14*1.5 for leaves
#r.setRmax_st([[12,9.0,6.0,12],[5.,5.],[15.]])
r.KMfu = 0.1#VERY IMPORTANT TO KEEP IT HIGH
#r.exud_k = np.array([2.4e-4])#*10#*(1e-1)
#r.k_gr = 1#0
r.sameVolume_meso_st = False
r.sameVolume_meso_seg = True
r.withInitVal =True
r.initValST = 0#0.6
r.initValMeso = 0#0.9
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


r.g0 = 8e-2
r.VcmaxrefChl1 =3
r.VcmaxrefChl2 = 10
r.a1 = 2#0.8/0.2
r.a3 = 2.2
r.Rd_ref = 2e-6
r.alpha = 0.27
r.theta = 0.51
krm2 = 2e-5 
krm1  =3e-03

rho_suc = {2 : 0.1,#0.8,#1.,#15,
           3 : 0.1,
           4 : 0.1}

maintenanceTot = 0
structSumInit = 0
orgs_all = r.plant.getOrgans(-1, True)
#print(len(orgs_all),len(r.plant.getOrgans(2, True)),len(r.plant.getOrgans(3, True)),len(r.plant.getOrgans(4, True)), r.plant.getNumberOfOrgans())
#print(len(orgs_all),len(r.plant.getOrgans(2)),len(r.plant.getOrgans(3)),len(r.plant.getOrgans(4,)), r.plant.getNumberOfOrgans())
#raise Exception
i__ = 0
for org in orgs_all:
    if org.organType() < 2:
        raise Exception("ot < 3")
    structSumInit += org.orgVolume(-1,False) * rho_suc[org.organType()]
    #print("num",i__, "vol",org.orgVolume(-1,False) ,"rho", rho_suc[org.organType()], "tot ",structSumInit)
    i__ +=1
AnSum = 0
i= 0

while simDuration < simMax: 
    
    print('simDuration:',simDuration )
    
    weatherX = weather(simDuration)

    r.Qlight = weatherX["Qlight"] #; TairC = weatherX["TairC"] ; text = "night"
        
        
    setKrKx_xylem(weatherX["TairC"], weatherX["RH"])
    
    #r.maxLoop = 3
    r.solve_photosynthesis(sim_time_ = simDuration, sxx_=sx, cells_ = True,RH_ = weatherX["RH"],
        verbose_ = False, doLog_ = False,TairC_= weatherX["TairC"] )
        
    AnSum += np.sum(r.Ag4Phloem)*dt
    
    startphloem= simDuration
    endphloem = startphloem + dt
    stepphloem = 1    
    errLeuning = sum(r.outputFlux)
    fluxes = np.array(r.outputFlux)
    segIdx = r.get_segments_index(4)
    idleafBlade = np.where(np.array(r.plant.leafBladeSurface)[segIdx] > 0)
    r.ci_adapt = False
    print("An",np.mean(np.array(r.An)[idleafBlade])*1e6, "mumol CO2 m-2 s-1")#An.append
    print("Rd",r.Rd*1e6, "mumol CO2 m-2 s-1")#An.append
    print("Vc",np.mean(np.array(r.Vc)[idleafBlade])*1e6, "mumol CO2 m-2 s-1")#Vc.append
    #print("Vcmax",np.mean(r.Vcmax)*1e6, "mumol CO2 m-2 s-1")#Vc.append
    print("Vj",np.mean(np.array(r.Vj)[idleafBlade])*1e6, "mumol CO2 m-2 s-1")#Vj.append
    #print("Vjmax",np.mean(r.Vjmax)*1e6, "mumol CO2 m-2 s-1")#Vj.append
    print("gco2",np.mean(np.array(r.gco2)[idleafBlade]), "mol CO2 m-2 s-1")#gco2.append
    print("gh2o",np.mean(np.array(r.gco2)[idleafBlade])*1.6, "mol H2O m-2 s-1")#gco2.append
    print("cics",np.mean(np.array(r.ci)[idleafBlade])/r.cs,"mol mol-1")#cics.append
    print("ci",np.mean(np.array(r.ci)[idleafBlade]), "mol mol-1")#fw.append
    print("deltagco2",np.mean(np.array(r.deltagco2)[idleafBlade]), "mol mol-1")#fw.append
    print("fw",np.mean(np.array(r.fw)[idleafBlade]), "-")#fw.append
    print("fluxes ", sum(fluxes))
    segRootIdx = r.get_segments_index(2)
    print("trans ",sum(fluxes[segIdx]),sum(fluxes[segRootIdx]), "cm3/day")
    #print(fluxes[segIdx],r.plant.leafBladeSurface)
    print("trans ",sum(fluxes[segIdx]),sum(r.plant.leafBladeSurface))
    #print(np.array(r.fw)[idleafBlade])
    #print(np.array(r.ci)[idleafBlade])
    #print(np.array(r.plant.leafBladeSurface)[segIdx][idleafBlade])
    ana = pb.SegmentAnalyser(r.plant.mappedSegments())
    
    cutoff = 1e-15 #is get value too small, makes paraview crash
    psiXyl_p = np.array(r.psiXyl)
    psiXyl_p[abs(psiXyl_p) < cutoff] = 0
    
    ana.addData("psi_Xyl",psiXyl_p)
    ana.write("results"+directoryN+"15pm_"+ str(ö) +".vtp", ["psi_Xyl","organType","subType"])
    
    ö +=1
    
    write_file_array("psiXyl", r.psiXyl)
    write_file_array("trans", r.Ev)
    write_file_array("transrate",r.Jw)
    
    verbose_simulate = False
    r.plant.simulate(dt, verbose_simulate)#, "outputpm15.txt") #time span in days /(60*60*24)
    r2.plant.simulate(dt,  verbose_simulate)#, "outputpm15.txt")
    
    orgs_all = r.plant.getOrgans(-1, True)
    structSum = 0
    for org in orgs_all: #cm3 * mmol Suc /cm3 = mmol Suc 
        structSum += org.orgVolume(-1,False) * rho_suc[org.organType()]
    GrowthSum = structSum - structSumInit    
    simDuration += dt
    maintenanceTot += structSum*krm1
    Y = 1
    print("GrowthSum", GrowthSum*(1/Y),"resp ",AnSum*0.95 - GrowthSum*(1/Y), maintenanceTot, (AnSum*0.95 - GrowthSum*(1/Y))/AnSum*100, "exud ", AnSum*0.05, "mmol Suc", "init suc",structSumInit, "end sum", structSum)

structSum = 0
orgs_all = r.plant.getOrgans(-1, True)

for org in orgs_all: #cm3 * mmol Suc /cm3 = mmol Suc 
    structSum += org.orgVolume(-1,False) * rho_suc[org.organType()]
    

print("simDuration", simDuration, "d")
GrowthSum = structSum - structSumInit
print("simDuration", simDuration, "d")
print("AnSum", AnSum, "mmol Suc")
print("GrowthSum", GrowthSum*(1/Y),"resp ",AnSum*0.95 - GrowthSum*(1/Y), maintenanceTot, (AnSum*0.95 - GrowthSum*(1/Y))/AnSum*100, "exud ", AnSum*0.05, "mmol Suc", "init suc",structSumInit, "end sum", structSum)

end = datetime.now()
print(end - beginning)