""" water movement within the root (static soil) """
import sys; sys.path.append("../.."); sys.path.append("../../src/python_modules")
from photosynthesis_cpp import PhotosynthesisPython  # Python hybrid solver
#from Leuning_speedup import Leuning #about 0.7 for both
#from photosynthesis_cython import Leuning
import plantbox as pb
#from plantbox import Photosynthesis as ph
import vtk_plot as vp
import math
import numpy as np

#import matplotlib.pyplot as plt
from datetime import datetime, timedelta

""" Parameters """
hPa2cm = 1.0197
TairC = 22
TairK = TairC + 273.15

#mg/cm3

def setKrKx_xylem(TairC): #inC
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
    beta = 0.9    
    kz_l  = VascBundle_leaf *(rad_x_l_1 + rad_x_l_2)    *np.pi /(mu * 8)  
    kz_s  = VascBundle_stem *(rad_x_s_1 + rad_x_s_2)    *np.pi /(mu * 8) 
    kz_r0 = VascBundle_root * rad_x_r0_1                *np.pi /(mu * 8)  
    kz_r1 = VascBundle_root *(rad_x_r12_1 + rad_x_r12_2)*np.pi /(mu * 8) 
    kz_r2 = VascBundle_root *(rad_x_r12_1 + rad_x_r12_2)*np.pi /(mu * 8)  
    kz_r3 = VascBundle_root * rad_x_r3_1                *np.pi /(mu * 8) # 4.32e-1

    #radial conductivity [1/day],
    kr_l  = 3.83e-4 * hPa2cm# init: 3.83e-4 cm/d/hPa
    kr_s  = 0#1.e-20  * hPa2cm # set to almost 0
    kr_r0 = 6.37e-5 * hPa2cm 
    kr_r1 = 7.9e-5  * hPa2cm 
    kr_r2 = 7.9e-5  * hPa2cm  
    kr_r3 = 6.8e-5  * hPa2cm 
    
    r.setKr([[kr_r0,kr_r1,kr_r2,kr_r3],[kr_s,kr_s ],[kr_l]], kr_length = 0.8) 
    r.setKx([[kz_r0,kz_r1,kz_r2,kz_r3],[kz_s,kz_s ],[kz_l]])
    
    
    Rgaz=8.314 #J K-1 mol-1 = cm^3*MPa/K/mol
    rho_h2o = dEauPure/1000#g/cm3
    Mh2o = 18.05 #g/mol
    MPa2hPa = 10000
    hPa2cm = 1/0.9806806
    #log(-) * (cm^3*MPa/K/mol) * (K) *(g/cm3)/ (g/mol) * (hPa/MPa) * (cm/hPa) =  cm                      
    p_a = np.log(RH) * Rgaz * rho_h2o * (TairC + 273.15)/Mh2o * MPa2hPa * hPa2cm

    r.air_psi = p_a #*MPa2hPa #used only with xylem


def checkBounds(r_):
    orgs_root = r_.plant.getOrgans(2)
    minZ = np.Inf
    maxZ = -np.Inf
    for org in orgs_root:
        numnodes = org.getNumberOfNodes()
        n_ = 0
        while n_ < numnodes:
            nd = org.getNode(n_)
            n_ += 1
            if nd.z < minZ:
                minZ = nd.z
            if nd.z > maxZ:
                maxZ = nd.z
            if nd.z > 0:
                print("root too high ",nd.z)
                raise Exception
    print("Deepest root depth (cm):", minZ)  
    print("Highest root depth (cm):", maxZ)
    orgs_stem = r_.plant.getOrgans(3)
    minZ = np.Inf
    maxZ = -np.Inf
    orgs_leaf = r_.plant.getOrgans(4)
    for org in orgs_leaf:
        numnodes = org.getNumberOfNodes()
        n_ = 0
        while n_ < numnodes:
            nd = org.getNode(n_)
            n_ += 1
            if nd.z < 0:
                print("leaf too low ", nd.z)
                #raise Exception
            if nd.z < minZ:
                minZ = nd.z
            if nd.z > maxZ:
                maxZ = nd.z
    print("Deepest root depth (cm):", minZ)  
    print("Highest root depth (cm):", maxZ)

# gmax = 0.004 #  cm3/day radial conductivity between xylem and guard cell
r, depth = 8/2, 60 # Soil core analysis
#p_s = -1000.  # static water potential (saturation) 33kPa in cm
p_mean = -187 #for theta = 0.35
p_bot = p_mean + depth/2
p_top = p_mean - depth/2
p_s = np.linspace(p_top, p_bot, depth)
#p_g = -2000 # water potential of the guard cell
RH = 0.45 # relative humidity
#TairC = 20
#p_a =  -10000  #default outer water potential 
k_soil = []
cs = 750e-6 #co2 paartial pressure at leaf surface (mol mol-1)
#TairK = TairC + 273.15


es = 0.61078 * math.exp(17.27 * TairC / (TairC + 237.3)) /10 #hPa
ea = es * RH 
VPD = es - ea #hPa

# root system 
pl = pb.MappedPlant(seednum = 2) #pb.MappedRootSystem() #pb.MappedPlant()
path = "../../modelparameter/plant/" #"../../../modelparameter/rootsystem/" 
name ="Triticum_aestivum_adapted_2021"#"wheat_uqr15" #"Triticum_aestivum_adapted_2021"#"wheat_uqr15" #"manyleaves" # "Anagallis_femina_Leitner_2010"  # Zea_mays_1_Leitner_2010
pl.readParameters(path + name + ".xml")

""" soil """

soilSpace =  pb.SDF_PlantContainer(r, r, depth, False)#
pl.setGeometry(soilSpace) # creates soil space to stop roots from growing out of the soil
pl.initialize(verbose = True)

# for rp11 in range(len(pl.organParam[2])):
    # print(rp11,"f_tf",pl.organParam[2][rp11].f_tf.alphaN )
    # if(pl.organParam[2][rp11].f_tf.alphaN <40):
        # pl.organParam[2][rp11].f_tf.alphaN +=20
        # print(rp11,"f_tf",pl.organParam[2][rp11].f_tf.alphaN )
        
simtime =7# [day] for task b
#Q = 770e-6#900e-6 # mol quanta m-2 s-1 light, example from leuning1995
Qs = np.array([0,580e-6]) # 580e-6, mol quanta m-2 s-1 light
pl.simulate(simtime)
#raise Exception
soil_index = lambda x,y,z : max(int(np.floor(-z)),-1) #abovegroud nodes get index -1
pl.setSoilGrid(soil_index)

r = PhotosynthesisPython(pl,psiXylInit = p_top,ciInit = cs*0.6)#PhloemFluxPython, pb.Photosynthesis(pl) 
#checkBounds(r)

SPAD= 41
print(0.0599*np.exp(0.0493*SPAD))
print((0.114 *(SPAD**2)+ 7.39 *SPAD+ 10.6)/10) #microg cm-2
chl_ = (0.114 *(SPAD**2)+ 7.39 *SPAD+ 10.6)/10
r.Chl = np.array( [chl_]) 
r.cs = cs
#r.g0 = 0.025
beginning = datetime.now()
organTypes = np.array(r.rs.organTypes)
trans = []
An = []
Vc = []
Vj = []
gco2 = []
cics = []
fw = []
pl = []
r.g0 = 8e-3
r.VcmaxrefChl1 =2#1.3#1.28#/2
r.VcmaxrefChl2 = 9#8.4#8.33#/2
r.a1 = 0.5#0#0.6/0.4#0.7/0.3#0.6/0.4 #ci/(cs - ci) for ci = 0.6*cs
r.a3 = 1.5
r.alpha = 0.4#0.2#/2
r.theta = 0.6#0.9#/2

# print(r.plant.dist2tips)
setKrKx_xylem(TairC)
# #Additional vtk plot
ana = pb.SegmentAnalyser(r.plant)
ana.write("results/photo28.vtp")#, ["radius", "surface", "rx", "fluxes"]) #
                                                                        
 #vp.plot_roots(ana, "rx", "Xylem matric potential (cm)")  # "fluxes"
 
trans28 = np.array([0.,0.])
trans28bis = np.array([0.,0.])
test = r.plant.getSegmentIds(4)
print(r.plant.getSegmentIds(-1))
print(len(r.plant.getSegmentIds(2)))
print(len(r.plant.getSegmentIds(3)))
print(len(r.plant.getSegmentIds(4)))

for i, Q_ in enumerate(Qs):
    print(Q_)
    r.Qlight = Q_
    r.maxLoop = 6
    #print(r.plant.segExchangesurf)
    r.solve_photosynthesis(sim_time_ = simtime,sxx_=p_s, cells_ = True,RH_ = RH,verbose_ = 1, doLog_ = True,TairC_=TairC)
    fluxes = np.array(r.outputFlux) # cm3/day
    
    trans.append(sum(np.where(organTypes == 4, fluxes,0)))
    #mmol Suc d-1 * (1/cm2 leaf) * (12 mol Co2/ mol Suc) * (1e3 mumol/mmol) * (1000 cm2/m2) * ( 1/(24*60*60) d/s)
    #=mmol Suc * (/cm2 leaf) * ( mol Co2/ mol Suc) * ( cm2/m2) * ( /s)
    #to mumol CO2 m-2 leaf s-1
    #print("An_bis",np.mean(np.array(r.Ag4Phloem)/np.array(segExchangesurf))*12*1e3*( 1/(24*60*60)),"mumol CO2 m-2 leaf s-1")
    an28 = np.mean(r.An)*1e6
    print("An",np.mean(r.An)*1e6)#An.append
    print("Vc",np.mean(r.Vc)*1e6)#Vc.append
    print("Vj",np.mean(r.Vj)*1e6)#Vj.append
    print("gco2",np.mean(r.gco2))#gco2.append
    print("cics",np.mean(r.ci)/cs)#cics.append
    print("fw",np.mean(r.fw))#fw.append
    print("fluxes ", sum(fluxes))
    segIdx = r.get_segments_index(4)
    segRootIdx = r.get_segments_index(2)
    fluxesSoil = r.soilFluxes( simtime, r.psiXyl,p_s, approx=False)
    #fluxescheck = r.
    print("trans ",sum(fluxes[segIdx]),sum(fluxes[segRootIdx]), sum(fluxesSoil.values()), sum(fluxes[segRootIdx])-sum(fluxesSoil.values()))# , fluxesSoil)
    print(sum(r.sumSegFluxes(fluxes).values()))
    trans28[i] = sum(fluxes[segIdx])
    trans28bis[i] = sum(fluxesSoil.values())
    print("trans ",sum(fluxes[segIdx]), sum(fluxesSoil.values()))

# #raise Exception
raise Exception
# Qs = np.array([0,1500e-6]) # 580e-6, mol quanta m-2 s-1 light

# p_mean = -295
# p_bot = p_mean + depth/2
# p_top = p_mean - depth/2
# p_s = np.linspace(p_top, p_bot, depth)

# SPAD= 41.5
# print(0.0599*np.exp(0.0493*SPAD))
# print((0.114 *(SPAD**2)+ 7.39 *SPAD+ 10.6)/10) #microg cm-2
# chl_ = (0.114 *(SPAD**2)+ 7.39 *SPAD+ 10.6)/10
# r.Chl = np.array( [chl_]) 
# r.plant.simulate(35-28)
# checkBounds(r)
# ana = pb.SegmentAnalyser(r.plant)
# ana.write("results/photo35.vtp")#, ["radius", "surface", "rx", "fluxes"]) #

# trans35 = np.array([0.,0.])
# trans35bis = np.array([0.,0.])
# for i, Q_ in enumerate(Qs):
    # r.Qlight = Q_
    # #print(r.plant.segExchangesurf)
    # r.solve_photosynthesis(sim_time_ = simtime,sxx_=p_s, cells_ = True,RH_ = RH,verbose_ = False, doLog_ = False,TairC_=TairC)
    # fluxes = np.array(r.outputFlux) # cm3/day
    
    # #trans.append(sum(np.where(organTypes == 4, fluxes,0)))
    # #mmol Suc d-1 * (1/cm2 leaf) * (12 mol Co2/ mol Suc) * (1e3 mumol/mmol) * (1000 cm2/m2) * ( 1/(24*60*60) d/s)
    # #=mmol Suc * (/cm2 leaf) * ( mol Co2/ mol Suc) * ( cm2/m2) * ( /s)
    # #to mumol CO2 m-2 leaf s-1
    # #print("An_bis",np.mean(np.array(r.Ag4Phloem)/np.array(segExchangesurf))*12*1e3*( 1/(24*60*60)),"mumol CO2 m-2 leaf s-1")
    # an35 = np.mean(r.An)*1e6
    # print("An",np.mean(r.An)*1e6)#An.append
    # print("Vc",np.mean(r.Vc)*1e6)#Vc.append
    # print("Vj",np.mean(r.Vj)*1e6)#Vj.append
    # print("gco2",np.mean(r.gco2))#gco2.append
    # print("cics",np.mean(r.ci)/cs)#cics.append
    # print("fw",np.mean(r.fw))#fw.append
    # print("fluxes ", sum(fluxes))
    # segIdx = r.get_segments_index(4)
    # fluxesSoil = r.soilFluxes( simtime, r.psiXyl,p_s, approx=False)
    # print("trans ",sum(fluxes[segIdx]), sum(fluxesSoil))
    # trans35[i] = sum(fluxes[segIdx])
    # trans35bis[i] = sum(fluxesSoil.values())


# Qs = np.array([0,1500e-6]) # 580e-6, mol quanta m-2 s-1 light
# p_mean = -263
# p_bot = p_mean + depth/2
# p_top = p_mean - depth/2
# p_s = np.linspace(p_top, p_bot, depth)
# SPAD= 38
# print(0.0599*np.exp(0.0493*SPAD))
# print((0.114 *(SPAD**2)+ 7.39 *SPAD+ 10.6)/10) #microg cm-2
# chl_ = (0.114 *(SPAD**2)+ 7.39 *SPAD+ 10.6)/10
# r.Chl = np.array( [chl_]) 
# r.plant.simulate(42-35)
# checkBounds(r)
# ana = pb.SegmentAnalyser(r.plant)
# ana.write("results/photo42.vtp")#, ["radius", "surface", "rx", "fluxes"]) #

# trans42 = np.array([0.,0.])
# trans42bis = np.array([0.,0.])
# for i, Q_ in enumerate(Qs):
    # print(Q_)
    # r.Qlight = Q_
    # #print(r.plant.segExchangesurf)
    # r.solve_photosynthesis(sim_time_ = simtime,sxx_=p_s, cells_ = True,RH_ = RH,verbose_ = 2, doLog_ = False,TairC_=TairC)
    # fluxes = np.array(r.outputFlux) # cm3/day
    
    # #trans.append(sum(np.where(organTypes == 4, fluxes,0)))
    # #mmol Suc d-1 * (1/cm2 leaf) * (12 mol Co2/ mol Suc) * (1e3 mumol/mmol) * (1000 cm2/m2) * ( 1/(24*60*60) d/s)
    # #=mmol Suc * (/cm2 leaf) * ( mol Co2/ mol Suc) * ( cm2/m2) * ( /s)
    # #to mumol CO2 m-2 leaf s-1
    # #print("An_bis",np.mean(np.array(r.Ag4Phloem)/np.array(segExchangesurf))*12*1e3*( 1/(24*60*60)),"mumol CO2 m-2 leaf s-1")
    # an42 = np.mean(r.An)*1e6
    # print("An",np.mean(r.An)*1e6)#An.append
    # print("Vc",np.mean(r.Vc)*1e6)#Vc.append
    # print("Vj",np.mean(r.Vj)*1e6)#Vj.append
    # print("gco2",np.mean(r.gco2))#gco2.append
    # print("cics",np.mean(r.ci)/cs)#cics.append
    # print("fw",np.mean(r.fw))#fw.append
    # print("fluxes ", sum(fluxes))
    # segIdx = r.get_segments_index(4)
    # #print("get fluxes soil", r.psiXyl[5074], r.psiXyl_old[5074], fluxes[5074], r.outputFlux_old[5074])
    # fluxesSoil = r.soilFluxes( simtime, r.psiXyl,p_s, approx=False)
    # #print("get fluxes soil_end")
    # print("trans ",sum(fluxes[segIdx]), sum(fluxesSoil))
    # trans42[i] = sum(fluxes[segIdx])
    # trans42bis[i] = sum(fluxesSoil.values())
    
    # # segRootIdx = r.get_segments_index(2)
    # # print("trans ",sum(fluxes[segIdx]),sum(fluxes[segRootIdx]), sum(fluxesSoil.values()), sum(fluxes[segRootIdx])-sum(fluxesSoil.values()))# , fluxesSoil)
    # # print(sum(r.sumSegFluxes(fluxes).values()))
    # # rx = r.psiXyl
    # # fluxes2 = np.array(r.segFluxes2(simtime,rx ,p_s,approx=False, cells= True))
    # # fluxes3 = np.array(r.segFluxes2(simtime, rx,p_s, approx=False, cells= True))
    # # fmin = np.minimum(np.absolute(fluxes2), np.absolute(fluxes))
    # # fmin[np.where(fmin == 0)]=1
    # # diif = np.absolute((fluxes2-fluxes)/fmin)
    # # print("diif ",max(diif), sum(diif), sum(np.array(fluxes2)- np.array(fluxes)))
    # # print("trans ",sum(fluxes[segIdx]),sum(fluxes2[segIdx]))
    # # print("trans ",sum(fluxes[segRootIdx]),sum(fluxes2[segRootIdx]))
    # # print("trans ",sum(fluxes[r.get_segments_index(3)]),sum(fluxes2[r.get_segments_index(3)]))
    # # print(sum(fluxes2), sum(fluxes3))
    # # r.solve_leuning_(sim_time_ = simtime,sxx_=p_s, cells_ = True,RH_ = RH,verbose_ = 1, doLog_ = False,TairC_=TairC)
    # # fluxes = np.array(r.outputFlux) # cm3/day
    # # print("trans ",sum(fluxes[segIdx]), sum(fluxes[segRootIdx]))

# # print("to check light saturation")
# # Qs = np.array([300e-6]) # 580e-6, mol quanta m-2 s-1 light
# # for Q_ in Qs:
    # # print(Q_)
    # # r.Qlight = Q_
    # # #print(r.plant.segExchangesurf)
    # # r.solve_leuning_(sim_time_ = simtime,sxx_=p_s, cells_ = True,RH_ = RH,verbose_ = False, doLog_ = False,TairC_=TairC)
    # # fluxes = np.array(r.outputFlux) # cm3/day
    
    # # #trans.append(sum(np.where(organTypes == 4, fluxes,0)))
    # # #mmol Suc d-1 * (1/cm2 leaf) * (12 mol Co2/ mol Suc) * (1e3 mumol/mmol) * (1000 cm2/m2) * ( 1/(24*60*60) d/s)
    # # #=mmol Suc * (/cm2 leaf) * ( mol Co2/ mol Suc) * ( cm2/m2) * ( /s)
    # # #to mumol CO2 m-2 leaf s-1
    # # #print("An_bis",np.mean(np.array(r.Ag4Phloem)/np.array(segExchangesurf))*12*1e3*( 1/(24*60*60)),"mumol CO2 m-2 leaf s-1")
    # # print("An",np.mean(r.An)*1e6)#An.append
    # # print("Vc",np.mean(r.Vc)*1e6)#Vc.append
    # # print("Vj",np.mean(r.Vj)*1e6)#Vj.append
    # # print("gco2",np.mean(r.gco2))#gco2.append
    # # print("cics",np.mean(r.ci)/cs)#cics.append
    # # print("fw",np.mean(r.fw))#fw.append
    # # print("fluxes ", sum(fluxes))
    # # segIdx = r.get_segments_index(4)
    # # print(sum(fluxes[segIdx]))
    # # #trans42 = sum(fluxes[segIdx])

# end = datetime.now()
# duration = end - beginning
# print('duration simulation ', duration, 's')

# print(trans28, trans35, trans42)
# print(trans28bis, trans35bis, trans42bis)
# print(an28, an35, an42)

