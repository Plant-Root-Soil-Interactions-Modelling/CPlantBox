""" water movement within the root (static soil) """

import os
main_dir=os.environ['PWD']#dir of the file
directoryN = "/"+os.path.basename(__file__)[:-3]+"/"
results_dir = main_dir +"/results"+directoryN

if not os.path.exists(results_dir):
    os.makedirs(results_dir)
else:
    test = os.listdir(results_dir)
    for item in test:
        try:
            os.remove(results_dir+item)
        except:
            pass

import sys; 
sys.path.append("../../.."); 
sys.path.append("../../../src/python_modules")
CPBdir = "../../.."

from phloem_flux import PhloemFluxPython  
import plantbox as pb
import math
import numpy as np
import vtk_plot as vp
from datetime import datetime, timedelta
import matplotlib
matplotlib.use('AGG') 
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
#from plotnine import *

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
    name2 = 'results'+ directoryN+ name+ '.txt'
    with open(name2, 'a') as log:
        log.write(','.join([num for num in map(str, data)])  +'\n')

def write_file_float(name, data):
    name2 = 'results' + directoryN+  name+ '.txt'
    with open(name2, 'a') as log:
        log.write(repr( data)  +'\n')

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
    r.setKr([[kr_r0,kr_r1,kr_r2,kr_r0,kz_r3,kz_r3],[kr_s,kr_s ],[kr_l]], kr_length_=l_kr) 
    r.setKx([[kz_r0,kz_r1,kz_r2,kz_r0,kz_r3,kz_r3],[kz_s,kz_s ],[kz_l]])
    
    
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
    
    r.setKr_st([[kr_r0,kr_r1 ,kr_r2 ,kr_r0,kr_r3,kr_r3],[kr_s,kr_s ],[kr_l]] , kr_length_= l_kr)
    r.setKx_st([[kz_r0,kz_r12,kz_r12,kz_r0,kr_r3,kr_r3],[kz_s,kz_s ],[kz_l]])
    
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
    r.setAcross_st([[Across_s_r0,Across_s_r12,Across_s_r12,Across_s_r0,Across_s_r3,Across_s_r3],[Across_s_s,Across_s_s],[Across_s_l]])
    #actually, don t use perimeter currently
    #r.setPerimeter_st([[0,Perimeter_s_r0,Perimeter_s_r12,Perimeter_s_r12,Perimeter_s_r0],[0,Perimeter_s_s,Perimeter_s_s],[0,Perimeter_s_l]])
    
""" Parameters """

weatherInit = weather(0)
simInit = 35
simDuration = simInit # [day] init simtime
simMax =simInit+1
depth = 60
dt = 1/24 #1h
verbose = True

# plant system 
pl = pb.MappedPlant(seednum = 2) #pb.MappedRootSystem() #pb.MappedPlant()
path = CPBdir+"/modelparameter/plant/"
#path = "../../../CPlantBox_test_files/params/"
name = 'P0_plant_C'
# name = "Triticum_aestivum_adapted_2021"#"wheat_uqr15" #"manyleaves"## "Anagallis_femina_Leitner_2010"  # Zea_mays_1_Leitner_2010

pl.readParameters(path + name + ".xml")

for p in pl.getOrganRandomParameter(pb.leaf):
    p.lb =  0 # length of leaf stem
    p.la,  p.lmax = 49.12433414, 49.12433414
    p.areaMax = 71.95670914  # cm2, area reached when length = lmax
    N = 100  # N is rather high for testing
    phi = np.array([-90,-80, -45, 0., 45, 90]) / 180. * np.pi
    l = np.array([49.12433414,1 ,1, 0.3, 1, 49.12433414]) #distance from leaf center
    p.tropismT = 0 #plagiotropism #1 #gravitropism # 6: Anti-gravitropism to gravitropism
    p.tropismN = 5
    p.tropismS = 0.05
    p.tropismAge = 5 #< age at which tropism switch occures, only used if p.tropismT = 6
p.createLeafRadialGeometry(phi, l, N)

for p in pl.getOrganRandomParameter(pb.stem):
    r= 1.128705967
    p.r = r
    p.lmax = (28-7)*r   
    p.lb =  5
#raise Exception
sdf = pb.SDF_PlantBox(np.Inf, np.Inf, depth )

pl.setGeometry(sdf) # creates soil space to stop roots from growing out of the soil


pl.initialize(verbose = True)#, stochastic = False)
pl.simulate(simDuration, False)#, "outputpm15.txt")
ot_ =np.array(pl.organTypes)
segments_ = np.array(pl.segLength())
vp.plot_plant(pl, 'type')

""" Coupling to soil """



min_b = [-3./2, -12./2, -61.]#distance between wheat plants
max_b = [3./2, 12./2, 0.]
cell_number = [6, 24, 61]#1cm3? 
layers = depth; soilvolume = (depth / layers) * 3 * 12
k_soil = []
initial = weatherInit["p_mean"]#mean matric potential [cm] pressure head

p_mean = -187
p_bot = p_mean + depth/2
p_top = p_mean - depth/2
sx = np.linspace(p_top, p_bot, depth)

#p_g = -2000 # water potential of the gua
picker = lambda x,y,z : max(int(np.floor(-z)),-1) #abovegroud nodes get index -1

#picker = lambda x, y, z: s.pick([x, y, z])    
pl.setSoilGrid(picker)  # maps segment


""" Parameters: photosynthesis """
#create object for photosynthesis and phloem
r = PhloemFluxPython(pl,psiXylInit = min(sx),ciInit = weatherInit["cs"]*0.5)

r.g0 = 8e-3 # minimal stomatal opening
r.VcmaxrefChl1 =1.28 #parameter for effect of N on assimilation
r.VcmaxrefChl2 = 8.33 #parameter for effect of N on assimilation
r.a1 = 0.5 #link An to g (kg1 in paper)
r.a3 = 1.5 #g_co2 to g_h2o (kg2 in paper)
r.alpha = 0.4   #influences Vj, alpha in paper
r.theta = 0.6   #influences Vj, omega in paper

r.cs = weatherInit["cs"] #external carbon concentration 
SPAD= 41.0
chl_ = (0.114 *(SPAD**2)+ 7.39 *SPAD+ 10.6)/10
r.Chl = np.array( [chl_]) #mean leaf chlorophyl content
""" Parameters: phloem """
setKrKx_phloem()

r.setKrm2([[2e-5]])
r.setKrm1([[10e-2]])#([[2.5e-2]])
r.setRhoSucrose([[0.51],[0.65],[0.56]])#0.51
r.setRmax_st([[14.4,9.0,6.0,14.4],[5.,5.],[15.]])
r.KMfu = 0.1
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

r.expression = 6
r.update_viscosity = True
r.solver = 1
r.atol = 1e-12
r.rtol = 1e-8

r.Csoil = 1e-4


""" for post processing """

ö=0
beginning = datetime.now()
AnSum = 0
results=[]
resultsAn=[]
resultsgco2=[]
resultsVc=[]
resultsVj=[]
resultscics=[]
resultsfw=[]
resultspl=[]

while simDuration < simMax: 
    
    print('simDuration:',simDuration )
    
    weatherX = weather(simDuration)

    r.Qlight = weatherX["Qlight"]
        
    #reset conductivity to water every time as depends on temperature    
    setKrKx_xylem(weatherX["TairC"], weatherX["RH"])
    
    #compute photosynthesis
    r.solve_photosynthesis(sim_time_ = simDuration, sxx_=sx, cells_ = True,RH_ = weatherX["RH"],
        verbose_ = False, doLog_ = False,TairC_= weatherX["TairC"] )
    
    #for post-processing: cummulative assimilated carbon
    AnSum += np.sum(r.Ag4Phloem)*dt
    organTypes = np.array(r.get_organ_types())#per node
    errLeuning = sum(r.outputFlux) #radial fluxes, should sum to 0
    fluxes = np.array(r.outputFlux)
    leavesSegs = np.where(organTypes==4)
    fluxes_leaves = fluxes[leavesSegs]
    if (min(r.Ev) < 0) or (min(r.Jw) < 0) or (min(fluxes_leaves)<0):
        print("leaf looses water", min(r.Ev),min(r.Jw), min(fluxes_leaves))
        raise Exception
    
    
    results.append(sum(np.where(organTypes == 4, fluxes,0))) #total leaf water exchange (transpiration)
    leafBlades = np.where(np.array(r.ci) > 0)[0]
    resultsAn.append(np.mean(np.array(r.An)[leafBlades])*1e6) #assimilation
    resultsVc.append(np.mean(np.array(r.Vc)[leafBlades])*1e6) #rate of carboxylation
    resultsVj.append(np.mean(np.array(r.Vj)[leafBlades])*1e6) #rate of photon flow
    resultsgco2.append(np.mean(np.array(r.gco2)[leafBlades])) # stomatal opening
    resultscics.append(np.mean(np.array(r.ci)[leafBlades])/r.cs) #ci/cs ratio
    resultsfw.append(np.mean(np.array(r.fw)[leafBlades])) #water scarcity factor for stomatal opening
    resultspl.append(np.mean(np.array(r.psiXyl)[leafBlades])) #leaf water potential
    
    """ soil water flow """   
    #in this example, we have a static soil
    fluxesSoil = r.soilFluxes(simDuration, r.psiXyl, sx, approx=False)
    #s.setSource(fluxesSoil.copy())  # richards.py 
    #s.solve(dt)
    #sx = s.getSolutionHead()  # richards.py    
    #min_sx, min_rx, max_sx, max_rx = np.min(sx), np.min(r.psiXyl), np.max(sx), np.max(r.psiXyl)
    #n = round((simDuration- simInit)/(simMax-simInit) * 100.)
    #print("[" + ''.join(["*"]) * n + ''.join([" "]) * (100 - n) + "], [{:g}, {:g}] cm soil [{:g}, {:g}] cm root at {:g} days {:g}"
     #       .format(min_sx, max_sx, min_rx, max_rx, s.simTime, r.psiXyl[0]))
          
    startphloem= simDuration
    endphloem = startphloem + dt
    stepphloem = 1
    filename = "results/pmincpb_" + str(simDuration) + "_15pm.txt" 
    r.startPM(startphloem, endphloem, stepphloem, ( weatherX["TairC"]  +273.15) , verbose_phloem, filename)
    
    Q_ST    = np.array(r.Q_out[0:Nt])
    Q_meso  = np.array(r.Q_out[Nt:(Nt*2)])
    Q_Rm    = np.array(r.Q_out[(Nt*2):(Nt*3)])
    Q_Exud  = np.array(r.Q_out[(Nt*3):(Nt*4)])
    Q_Gr    = np.array(r.Q_out[(Nt*4):(Nt*5)])  
    
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
    
    
    Q_Rmmax       = np.array(r.Q_out[(Nt*5):(Nt*6)])
    Q_Grmax       = np.array(r.Q_out[(Nt*6):(Nt*7)])
    Q_Exudmax     = np.array(r.Q_out[(Nt*7):(Nt*8)])
    
    Q_ST_i        = Q_ST      - Q_STbu
    Q_Par_i       = Q_out     - Q_Parbu
    Q_Rm_i        = Q_Rm      - Q_Rmbu
    Q_Gr_i        = Q_Gr      - Q_Grbu
    
    
    Q_Exud_i      = Q_Exud    - Q_Exudbu
    Q_meso_i      = Q_meso    - Q_mesobu
    
    Q_Rmmax_i     = Q_Rmmax   - Q_Rmmaxbu
    Q_Grmax_i     = Q_Grmax   - Q_Grmaxbu
    Q_Exudmax_i   = Q_Exudmax - Q_Exudmaxbu
    
    Q_out_i       = Q_Rm_i    + Q_Exud_i      + Q_Gr_i
    Q_outmax_i    = Q_Rmmax_i + Q_Exudmax_i   + Q_Grmax_i
    """ for post processing """
    
    
    if verbose :
        print("\n\n\n\t\tat ", int(np.floor(simDuration)),"d", int((simDuration%1)*24),"h",  round(r.Qlight *1e6),"mumol m-2 s-1")
        print("Error in photos:\n\tabs (cm3/day) {:5.2e}".format(errLeuning))
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
        
    """ paraview output """    
    ana = pb.SegmentAnalyser(r.plant.mappedSegments())
    
    cutoff = 1e-15 #is get value too small, makes paraview crash
    fluxes_p = fluxes
    fluxes_p[abs(fluxes_p) < cutoff] = 0
    
    psiXyl_p = np.array(r.psiXyl)
    psiXyl_p[abs(psiXyl_p) < cutoff] = 0
    ana.addData("fluxes", fluxes_p)
    ana.addData("psi_Xyl",psiXyl_p)
    ana.write("results"+directoryN+"photo_"+ str(ö) +".vtp", ["fluxes","psi_Xyl"]) 
    
      
    ö +=1
    
    
    """ print to files """    
    write_file_array("fluxes", fluxes)
    write_file_array("psiXyl", r.psiXyl)
    
    verbose_simulate = False
    r.plant.simulate(dt, verbose_simulate)#, "outputpm15.txt") #time span in days /(60*60*24)
    
        
    
    simDuration += dt
