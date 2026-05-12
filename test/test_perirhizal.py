
# this is a test file by Erik Kopp to test the alternative implementation of the perirhizal resistances and the perirhizal diffusion


import sys; sys.path.append("../src/functional")
sys.path.append("../../dumux-rosi/build-cmake/cpp/python_binding/")
sys.path.append("../../../dumuxtest/dumux/dumux-rosi/build-cmake/cpp/python_binding/")
#from plantbox import Perirhizal
from plantbox.functional.Perirhizal import PerirhizalPython
from numpy import linalg as LA
import plantbox.functional.van_genuchten as vg
from rosi.richards_flat import RichardsFlatWrapper as RichardsWrapper  # Python part, macroscopic soil model
#from richards import RichardsWrapper  # Python part
from rosi.richards_no_mpi import RichardsNoMPIWrapper  # Python part of cylindrcial
from rosi.rosi_richardsnc_cyl import RichardsNCCylFoam as RichardsNC_cyl # C++ part (Dumux binding), macroscopic soil model
#from rosi_richardsnc_cyl import RichardsNCCylFoam as RichardsNC_cyl  # C++ part (Dumux binding)
from rosi_richards22c import RichardsNCSPILU as RichardsNCSP #test

import Perirhizal
import pandas as pd
import numpy as np
import time

# run the dumux implementation of root water and nitrate uptake an then compare it to the alpha omega model

n_tests = 10 #try everything here for this many random parameter sets
do_computation = True #should the computation be run or take the data from a saved file

# general parameters

max_time = 10 #d
n_times = 100 # number of time intervals
times = np.linspace(0,max_time,n_times)[1:]
r_prhiz = 1 # cm
r_root = 0.02 # cm
rho = r_prhiz / r_root
NC = 101 # number of spatial discretisations
n_sp = NC - 1

dt = max_time / n_times

#space for the oupputs
watercontent_dumux = np.zeros((n_tests, n_times, n_sp+1))
soluteconcentration_dumux = np.zeros((n_tests, n_times, n_sp+1))
soluteconcentration_steadystate = np.zeros((n_tests, n_times, n_sp+1))
watercontent_steadyrate = np.zeros((n_tests, n_times, n_sp+1))
soluteconcentration_steadyrate = np.zeros((n_tests, n_times, n_sp+1))

soilVG = [0.078, 0.43, 0.036, 1.56, 24.96]  # hydrus loam soil

def load_constants(cyl):

    cyl.molarMassC = 12.011
    cyl.molarMassN = 14.0
    cyl.yr_per_d = 1/365 # [yr/d]
    cyl.m3_per_cm3 = 1e-6; # m3/cm3
    cyl.cm3_per_m3 = 1e6; # cm3/m3
    cyl.cm3_per_L = 1e3; # cm3/dm3

    cyl.kgC_per_mol = (1/1000) * cyl.molarMassC
    cyl.kgN_per_mol = (1/1000) * cyl.molarMassN
    
    yr_per_d = 1/365 # [yr/d]
    m3_per_cm3 = 1e-6; # m3/cm3
    cm3_per_m3 = 1e6; # cm3/m3

    molMarssW = 18 # g/mol
    rhoW = 1 #g/cm3
    rhoWM = rhoW/molMarssW #[g/cm3]*[g/mol] = mol/cm3
    cyl.rhoWM = rhoWM
    mlFr=1/rhoWM

    cyl.molarMassC = 12.011
    cyl.molarMassN = 14.0
    cyl.mg_per_molC = cyl.molarMassC * 1000.
    cyl.mg_per_molN = cyl.molarMassN * 1000.
    yr_per_d = 1/365 # [yr/d]
    m3_per_cm3 = 1e-6; # m3/cm3
    cm3_per_m3 = 1e6; # cm3/m3

    molMarssW = 18 # g/mol
    molMassMulC = 30 #g/mol
    rhoW = 1 #g/cm3
    rhoWM = rhoW/molMarssW #[g/cm3]*[g/mol] = mol/cm3
    cyl.mlFr=1/rhoWM #only water phase is relevant right now




    return cyl

def getSoilTextureAndShape(s, soilVG):
    """ soil shape and texture data
        to adapt according to the soil represented
    """
    
    constants = {'solidDensity': 2650, # kg/m3
    'solidMolarMass': 60.08e-3 #kg/mol 
    }
    paramSet = pd.DataFrame()
    
    solidDensity = 2650 # [kg/m^3 solid] #taken from google docs TraiRhizo
    solidMolarMass = 60.08e-3 # [kg/mol]
    # theta_r, theta_s, alpha, n, Ks
    #soilVG = [paramSet['soilVG_theta_r'],paramSet['soilVG_theta_s'],paramSet['soilVG_alpha'],paramSet['soilVG_n'],paramSet['soilVG_Ks']]
    soilTextureAndShape = {
                           "solidDensity":solidDensity,
                        'solidMolarMass': solidMolarMass,
                           'soilVG':soilVG}
    #common transformations
    
    soilTexture =soilTextureAndShape
    
    s.solidDensity = soilTexture['solidDensity'] #[kg/m^3 solid]
    s.solidMolarMass = soilTexture['solidMolarMass']# [kg/mol]
    s.soil =  soilTexture['soilVG']
    
    s.vg_soil = vg.Parameters(s.soil)
    # [mol / m3 solid] =[kg/m^3 solid] / [kg/mol]
    s.solidMolDensity = s.solidDensity/s.solidMolarMass
    # [mol / m3 scv] = [mol / m3 solid] * [m3 solid /m3 space]
    s.bulkDensity_m3 = s.solidMolDensity*(1.- s.vg_soil.theta_S)
    # (kg/m3) * (volume soil / total volume) * (g/kg) * (m3/cm3)
    s.bulkMassDensity_gpercm3 = s.solidDensity*(1.- s.vg_soil.theta_S)*1000*1e-6
    
    s.masspercm3 = s.solidDensity * 1.e-6
    return soilTexture

def run_perirhizal_test():
    
    #space for the outputs
    water_dumux = np.zeros((n_times, n_sp+1))
    solute_dumux = np.zeros((n_times, n_sp+1))
    water_ss = np.zeros((n_times, n_sp+1))
    solute_ss = np.zeros((n_times, n_sp+1))
    water_sr = np.zeros((n_times, n_sp+1))
    solute_sr = np.zeros((n_times, n_sp+1))
    
    mean_water = np.zeros((n_times))
    mean_solutes = np.zeros((n_times))
    
    # determine some random parameters
    initial_waterpotential = -100 + np.random.rand() * 50 #cm3/cm3 #or choose an initial pressure head?
    initial_soluteconcentration = 1e-6*(1.0+np.random.rand()) #mol/cm3 #TODO: lookup realistic concnetration, maybe 10 times the Michaelis Menten half saturation?
    print("initial_soluteconcentration",initial_soluteconcentration)
    
    # root conductivity and solute uptake parameters, it is chosen to be constant throughout the entire simulation time
    root_conductivity = 1e-4 #1/d
    inner_kr = root_conductivity * r_root
    Vmax = 1.0e-10 #mol/d
    Km = 1.0e-9 #mol/cm3
    
    DS_W = 1.902e-5 #cm2/s
    Ds = DS_W * 3600 * 24 #cm2/d

    # the xylem matrix potential varies over time (keep it low so that there is little to no outflow of water)
    rx_t = lambda t : -500+200*np.sin(t) #cm

    # load the perirhizal model
    peri = PerirhizalPython() 
    sp = vg.Parameters(soilVG)
    peri.set_soil(sp)
    #no lookup tables are used here as there arent many simulations

    lb = 0.5
    points = np.logspace(np.log(r_root) / np.log(lb), np.log(r_prhiz) / np.log(lb),
                                  NC, base = lb)
    CC = np.array([(points[i] + points[i+1])/2 for i in range(len(points)-1)])

    

    # initialise the dumux model
    print("test0")
    #cyl = RichardsNoMPIWrapper(RichardsNCCylFoam())
    cyl = RichardsWrapper(RichardsNC_cyl())
    cyl = load_constants(cyl)
    soilTexture = getSoilTextureAndShape(cyl, soilVG)
    print("test005")
   # cyl = RichardsNoMPIWrapper(RichardsNC_cyl())
    c2 = [cyl.mlFr * CC[i] for i in range(NC-1)]
    #cyl = RichardsWrapper(RichardsCylFoam())
    #cyl.initialize()
    cyl.createGrid1d(points)
    initial=-10
    cyl.setHomogeneousIC(initial)  # cm pressure head
    #cyl.setICZ_solute(cyl.dumux_str([q for q in c2]))
    cyl.setICZ_solute(str([q for q in c2]))
    print("test006")
    cyl.initialize()
    
    cyl.setVGParameters([soilVG])
    cyl.setOuterBC("fluxCyl", 0.)  # [cm/day]
    cyl.setInnerBC("fluxCyl", 0.)  # [cm/day]
    print("test007")
    #cyl.initializeProblem()
    cyl.maxDt = 3500/(3600*24) # soilModel.maxDt_1DS
    #cyl.setParameter("Soil.mucilCAffectsW", "true")
    cyl.setParameter("Problem.verbose", "0")
    cyl.setParameter("Newton.Verbosity", "0")
    cyl.initializeProblem(maxDt=cyl.maxDt ) 
    print("test008")
    critP=-15000
    cyl.setCriticalPressure(critP)  # cm pressure head
    #cyls.append(cyl)
    print("test0.1")
    seg_length = 1 #cm
    
    if True:
        cyl.setParameter("Newton.EnableResidualCriterion", "true") # sometimes helps, sometimes makes things worse
        cyl.setParameter("Newton.EnableAbsoluteResidualCriterion", "true")
        cyl.setParameter("Newton.SatisfyResidualAndShiftCriterion", "true")
        cyl.setParameter("Newton.MaxRelativeShift", "1e-9")# reset value
    #cyl.initialize(verbose = False) # No parameter file found. Continuing without parameter file.
    #cyl.initialize()
    print("test1")
    #cyl.initializeProblem(maxDt=cyl.maxDt )
    cyl = load_constants(cyl)
    soilTexture = getSoilTextureAndShape(cyl, soilVG)
    
    c2 = [cyl.mlFr * CC[i] for i in range(NC-1)]
    
    cyl.createGrid1d(points)# cm
    cyl.doSoluteFlow = True
    cyl.noAds = False
    
    #set soil parameters
    cyl.solidDensity = soilTexture['solidDensity'] #[kg/m^3 solid]
    cyl.solidMolarMass = soilTexture['solidMolarMass']# [kg/mol]
    cyl.soil =  soilTexture['soilVG']

    cyl.setParameter( "Soil.MolarMass", str(cyl.solidMolarMass))
    cyl.setParameter( "Soil.solidDensity", str(cyl.solidDensity))
    #cyl.setVGParameters([cyl.soil])

    minD=1.0e-12#minimal diffusion coefficient, I think this will be numerically good
    #cyl.setParameter("1.Component.LiquidDiffusionCoefficient", str(DS_W/10000+minD)) #cm^2/s -> m^2/s
    #todo: reactivate the diffusion coefficient
    
    #set michaelis menten uptake parameters
    cyl.setParameter("RootSystem.Uptake.Vmax", str(Vmax))
    cyl.setParameter("RootSystem.Uptake.Km", str(Km))
    
    #cyl = setSoilParam(cyl,paramSet) # here the soil parameters are set out of paramSet
    #cyl = setBiochemParam(cyl,paramSet) # here the biochemical parameters are set out of paramSet
    print("test2") 
    #cyl.initializeProblem(maxDt=cyl.maxDt )    
    doCells = False
    if doCells:
        Cells = (cyl.getCellCenters_().reshape(-1))
    else:
        Cells = []
    cyl.setParameter("Flux.UpwindWeight", "1")
    cyl.setParameter("SpatialParams.Temperature","293.15") # todo: redefine at each time step
    cyl.setParameter("Soil.BC.dzScaling", "1")
    #cyl.setParameter( "Soil.css1Function", str(9)) # todo: can I activate this?
    cyl.setParameter("Problem.verbose", "0")
    cyl.setParameter("Problem.reactionExclusive", "0")
    cyl.setParameter("Soil.CriticalPressure", str(-15000))
    cyl.seg_length = seg_length
    cyl.setParameter("Problem.segLength", str(seg_length))   # cm
    cyl.l = seg_length
    cyl.setParameter( "Soil.Grid.Cells",str( NC-1)) # -1 to go from vertices to cell (dof) # I am unsure what this comment means
    cyl.setParameter("Newton.MaxTimeStepDivisions",
                str( 100) )
    cyl.setParameter("Newton.MaxSteps",
                str( 100) )
    cyl.setParameter("Newton.EnableChop", "true")
    print("test3")
    #cyl.initializeProblem(maxDt=cyl.maxDt )
    #boundary conditions
    #water
    cyl.setParameter( "Soil.BC.Bot.Type", str(int(3)))
    cyl.setParameter( "Soil.BC.Top.Type", str(int(3)))
    cyl.setParameter( "Soil.BC.Bot.Value", str(0.0)) #will be prescribed at each timestep
    cyl.setParameter( "Soil.BC.Top.Value", str(0.0))
    print("test3.1")
    #solutes
    cyl.setParameter( "Soil.BC.Bot.SType", str(int(3))) # TODO change this to 8 as Michaelis mnetne
    cyl.setParameter( "Soil.BC.Top.SType", str(int(3)))
    cyl.setParameter( "Soil.BC.Bot.CValue", str(0.0)) #Michaelis Menten parameters are set elsewhere
    cyl.setParameter( "Soil.BC.Top.CValue", str(0.0))
    
    print("test3.2")
    #cyl.initializeProblem(maxDt=cyl.maxDt )
    
    #initial conditions
    cyl.setHomogeneousIC(initial_waterpotential)
    #cyl.setParameter("Soil.IC.P", cyl.dumux_str(initial_waterpotential))# cm #this is where the initial matrix potential and with it the water content are set
    #molar_fraction = initial_soluteconcentration / 62 * 18 /1 # / molar mass NO3 * molarMass Water
    #print("molar_fraction",molar_fraction)
    #cyl.setICZ_solute(c=[molar_fraction])
    theta = 0.3 #change this after testing
    c2 = [cyl.mlFr * initial_soluteconcentration/cyl.mg_per_molC*cyl.masspercm3/theta for i in range(NC-1)]
    #mlFr*paramSet['Soil.IC.C1']/cyl.mg_per_molC*masspercm3/theta
    temp = np.array(c2)*initial_soluteconcentration
    #temp2 = cyl.dumux_str(temp.tolist())
    temp2 = cyl.dumux_str(c2)
    #cyl.initializeProblem(maxDt=cyl.maxDt )
    print("IC",temp2)
    #cyl.initialize()
    #cyl.initializeProblem(maxDt=cyl.maxDt )
    #
    #cyl.initialize()
    #cyl.setParameter("Soil.IC.C1",cyl.dumux_str([q for q in c2]))
    #cyl.setParameter("Soil.IC.C", temp2)# cm #does this work without spatial resolution
    #cyl.setParameter("Soil.IC.C", initial_soluteconcentration)# cm #does this work without spatial resolution
    
    print("test3.3")
    cyl.initializeProblem(maxDt=cyl.maxDt )
    
    if len(Cells) > 0:#in case we update a cylinder, to allow dumux to do interpolation
        assert(len(c2)==len(Cells))
        CellsStr = cyl.dumux_str(Cells/100)#cm -> m #changed by Erik
        cyl.setParameter("Soil.IC.Z",CellsStr)# m
        if len(Cells)!= len(x):
            print("Cells, x",Cells, x, len(Cells), len(x))
            raise Exception
         
        #j = 2
        for j in range( 1, cyl.numComp):
            cyl.setParameter("Soil.IC.C"+str(j)+"Z",CellsStr) # m
        if len(Cells)!= len( c2):
            print("Cells,  cAll[j-1]",Cells,  c2,
                    len(Cells), len(c2), j)
            raise Exception
    print("test3.4")
    cyl.maxDt = 3500/(3600*24) # soilModel.maxDt_1DS
    #cyl.setParameter("Soil.mucilCAffectsW", "true")
    cyl.setParameter("Problem.verbose", "0")
    cyl.setParameter("Newton.Verbosity", "0")
    cyl.initializeProblem(maxDt=cyl.maxDt ) 
        
    #cyl.setParameter("Problem.doNeffect","true") 

    cyl.eps_regularization = 1e-14
    #cyl.setRegularisation(cyl.eps_regularization, cyl.eps_regularization)
    #which one is the relative error?
    #pcEps    capillary pressure regularisation, krEps 	relative permeabiltiy regularisation
    #default 1e-6, 1e-6
    cyl.setRegularisation(cyl.eps_regularization, cyl.eps_regularization)

    cyl.setCriticalPressure(-15000)  # cm pressure head
    #cyl = setInitialConditions(cyl,paramSet)

    # run the dumux model
    #get water and solute next to the root
    print("test3.5")
    current_rs_potential = cyl.getSolutionHead() #todo rmeove
    print(current_rs_potential[0])
    print("test3.6")
    #current_rs_potential = current_rs_potential[0]
    current_rs_potential = initial_waterpotential
    print("test3.7")
    current_rs_concentration = cyl.getSolution(1) * cyl.rhoWM # todo remove
    print(current_rs_concentration)
    print(current_rs_concentration[0])
    print("test3.8")
    #current_rs_concentration = current_rs_concentration[0]
    current_rs_concentration = initial_soluteconcentration
    print("test4")
    for i in range(n_times):
        rx = rx_t(i*dt)
        #model root water uptake 
        root_wateruptake = root_conductivity * (rx - current_rs_potential)
        water_dumux[i,0] = root_wateruptake
        cyl.setParameter( "Soil.BC.Bot.Value", str(root_wateruptake))
        
        #model root solute uptake 
        root_soluteuptake = Vmax * current_rs_concentration * (Km + current_rs_concentration)
        solute_dumux[i,0] = root_wateruptake
        cyl.setParameter( "Soil.BC.Bot.Value", str(root_wateruptake)) #todo set BC solute
        print("test5")
        cyl.solve(dt)
        
        # save outputs of dumux
        print("test6")
        water_dumux[i,1:] = cyl.getWaterContent()
        solute_dumux[i,1:] = cyl.getSolution(0)
        
        mean_water = np.mean(water_dumux[i,1:])
        mean_solutes = np.mean(solutes_dumux[i,1:])
        
        
        current_rs_potential = cyl.getSolutionHead()
        current_rs_potential = current_rs_potential[0]
        current_rs_concentration = solute_dumux[i,1]

    

    
        # run the alpah omega model on waterflow (steady rate)
        # mean water content of the rhizosphere
        total_water = 0
        total_solute = 0
        for j in range(NC-1):
            total_water = total_water + water_dumux[i,j+1]*(points[j+1]**2 - points[j]**2)
            total_solute = total_solute + water_dumux[i,j+1]*solute_dumux[i,j+1]*(points[j+1]**2 - points[j]**2)
        mean_water = total_water / (points[NC-1]**2 - points[0]**2)
        total_solute = total_solute / (points[NC-1]**2 - points[0]**2)
        
        #translate the mean water content to a mean matrix potential
        sx = vg.pressure_head(mean_water, peri.sp)
        
        #start the steady rate solver
        h_sr = peri.soil_root_interface_(rx, sx, inner_kr, rho, peri.sp)
        
        #compute both matrix flux potentials:
        Phi_root = vg.fast_mfp[sp](h_sr)
        Phi_soil = vg.fast_mfp[sp](s_x)
        
        #compute the spatial watercontents
        Phi_A, Phi_C = peri.determine_mfp_function(Phi_root, Phi_soil, rho)
        #outer Phi
        Phi_out = Phi_A+Phi_C
        
        for j in range(NC-1):
            s = CC[j] / r_prhiz
            Phi = Phi_A*(s**2-np.log(s**2))+Phi_C
            water_sr[i,j] = vg.water_content(vg.fast_imfp[sp](),peri.sp)
        

        # run the alpha omega model on solute flow (both steady state and steady rate)
        waterflow = root_wateruptake
        result_solutes_ss = peri.soil_root_solutes_ss_(Phi_out, Phi_soil, total_solute, Vmax, Km, Ds, waterflow, peri.sp)
        result_solutes_sr = peri.soil_root_solutes_sr_(Phi_root, Phi_soil, rho, total_solute, Vmax, Km, Ds, waterflow, peri.sp)
        
        solute_ss[i,1] = result_solutes_ss
        solute_ss[i,0] = Vmax * result_solutes_ss / (Km + result_solutes_ss)
        solute_sr[i,1] = result_solutes_sr
        solute_sr[i,0] = Vmax * result_solutes_sr / (Km + result_solutes_sr)
        
        F0 = self.lookup_table_solutes((Phi_root[i],0)) # for the ratio of concentration next to the root to somewhere in the perirhizal zone
        D_tilde = 1/Ds/math.pow(sp.theta_S-sp.theta_R,13/3)*(sp.theta_S*sp.theta_S)
        for j in range(NC-2):
            s = CC[j+1] / r_prhiz
            Phi_current = Phi_A*(s**2-np.log(s**2))+Phi_C
            F = self.lookup_table_solutes((Phi_soil[i],0))-F0
            F_tilde=math.exp(D_tilde*F)
            solute_ss[i,1] = result_solutes_ss * F_tilde + (1-F_tilde) * solute_ss[i,0] / waterflow
            solute_sr[i,1] = result_solutes_sr * F_tilde + (1-F_tilde) * solute_sr[i,0] / waterflow
    

    return water_dumux, solute_dumux, solute_ss, water_sr, solute_sr

if do_computation:
    # save everything in the np arrays
    for i in range(n_tests):
        water_dumux, solute_dumux, solute_ss, water_sr, solute_sr = run_perirhizal_test()
        watercontent_dumux[i,:,:]=water_dumux[:,:]
        soluteconcentration_dumux[i,:,:]=solute_dumux[:,:]
        soluteconcentration_steadystate[i,:,:]=solute_ss[:,:]
        watercontent_steadyrate[i,:,:]=water_sr[:,:]
        soluteconcentration_steadyrate[i,:,:]=solute_sr[:,:]
    
    np.savez("test_perirhizal.npz", watercontent_dumux=watercontent_dumux, soluteconcentration_dumux=soluteconcentration_dumux, soluteconcentration_steadystate=soluteconcentration_steadystate, watercontent_steadyrate=watercontent_steadyrate, soluteconcentration_steadyrate=soluteconcentration_steadyrate)

else:
    simulation_results = np.load("test_perirhizal.npz")
    watercontent_dumux = simulation_results["watercontent_dumux"]
    soluteconcentration_dumux = simulation_results["soluteconcentration_dumux"]
    soluteconcentration_steadystate = simulation_results["soluteconcentration_steadystate"]
    watercontent_steadyrate = simulation_results["watercontent_steadyrate"]
    soluteconcentration_steadyrate = simulation_results["soluteconcentration_steadyrate"]


# compare both for the differint means of water / solute content
run = 5
timestep = 10

#plot water and nitrogen in that one voxel
fig, ax1 = figure_style.subplots12(1, 1)
ax2 = ax1.twinx()

water_dumux = watercontent_dumux[run, timestep, :]
water_perirhizal = watercontent_steadyrate[run, timestep, :]
solute_dumux = soluteconcentration_dumux[run, timestep, :]
solute_steadystate = soluteconcentration_steadystate[run, timestep, :]
solute_steadyrate = soluteconcentration_steadyrate[run, timestep, :]
linestyle_dumux = "solid"
linestyle_steadystate = "dotted"
linestyle_steadyrate = "dashed"
    
ax1.plot(times, water_dumux, "b", linestyle = linestyle_dumux, label = "water_dumux")
ax1.plot(times, water_perirhizal, "b", linestyle = linestyle_steadyrate, label = "water_perirhizal")
ax2.plot(times, solute_dumux, "m", linestyle = linestyle_dumux, label = "solute_dumux")
ax2.plot(times, solute_steadystate, "m", linestyle = linestyle_steadystate, label = "solute_steadystate")
ax2.plot(times, solute_steadyrate, "m", linestyle = linestyle_steadyrate, label = "solute_steadyrate")

ax1.set_xlabel("time (day)")
ax1.set_ylabel("water")
ax2.set_ylabel("nitrogen concentration")
ax1.legend(["watercontent cm3/cm3"], loc="upper left")
ax2.legend(["nitrogen concentration mol/cm3"], loc="upper right")

ax1.legend(loc="upper left")
ax2.legend(loc="upper right")
#np.save("input/" + filename + "_fp", np.vstack((sim_times_, -np.array(t_act_), np.array(q_soil_)))) 
plt.show()




# #---old stuff below


# # generate random inputs for one perirhizal segment

# """
        # rx             xylem matric potential [cm]
        # sx             bulk soil matric potential [cm]
        # inner_kr       root radius times hydraulic conductivity [cm/day] 
        # rho            geometry factor [1] (outer_radius / inner_radius)
        # c_sol          concentration of solutes in the bulk soil
        # Vmax           Michaelis Menten Kinetics maximal solute uptake rate
        # Km             Michaelis Menten Kinetics half saturation
        # sp             soil parameter: van Genuchten parameter set (type vg.Parameters)
        # sp_theta_r
        # sp_theta_s
        # sp_n
        # sp_alpha
        # sp_Ksat
# """
# labels = ["rx","sx","inner_kr","rho","c_sol","Vmax","Km","hsr","hsr_lookup","hsr_global","waterflow","sol_c","sol_c_ss","sol_c_sr","sol_U","sol_U_ss","sol_U_sr"]

# inputs = ["rx","sx","inner_kr","rho","c_sol","Vmax","Km"]
# outputs = ["hsr","hsr_lookup","hsr_global","sol_c","sol_c_ss","sol_c_sr","sol_U","sol_U_ss","sol_U_sr"]
# outputs = ["rx","sx","waterflow","sol_c","sol_c_ss","sol_c_sr","sol_U","sol_U_ss","sol_U_sr"]

# #kr is assumed to lie between .5 and 5e-6 cm/hPa/d = 0.5 to 5e-6 1/d 
# #r is assumed to be 0.03
# #multiply this number by 10 in order to not fall under the threshhold at which the soil hydraulic potential is chosen as a default

# Intervals = {
        # "rx": [-14000,-10],
        # "sx": [-1000,-10],
        # "inner_kr": [1.5e-7,1.5e-6], 
        # "rho": [10,199.0],
        # "c_sol": [1.0e-9,1.0e-5],
        # "Vmax": [1.0e-11,1.0e-10],
        # "Km": [1.0e-8,1.0e-6],
        # }
# Intervals = pd.DataFrame(Intervals, index = ["min","max"])

# ntests = 10 
# tests = pd.DataFrame(index = range(ntests),columns = labels)

# for one_label in Intervals.columns:
    # min_val = Intervals.loc["min",one_label]
    # max_val = Intervals.loc["max",one_label]

    # testvalues = (max_val - min_val) * np.random.rand(ntests) + min_val * np.ones(ntests)
    # tests[one_label] = testvalues

# peri = PerirhizalPython()
# #peri = PerirhizalPython(Perirhizal)


# rx = np.array(tests.loc[:,"rx"]) #* (-1)
# sx = np.array(tests.loc[:,"sx"]) #* (-1)
# inner_kr = np.array(tests.loc[:,"inner_kr"])
# rho = np.array(tests.loc[:,"rho"])
# c_sol = np.array(tests.loc[:,"c_sol"])
# Vmax = np.array(tests.loc[:,"Vmax"])
# Km = np.array(tests.loc[:,"Km"])
# Ds = 1.0e1
# r_root = 0.03

# hydrus_loam = [0.078, 0.43, 0.036, 1.56, 24.96]
# sp = vg.Parameters(hydrus_loam)
# peri.set_soil(sp)

# #create the lookup table
# #t = time.time()
# #peri.create_lookup("results/hydrus_loam", sp)
# #elapsed_time = time.time() - t
# #print("Creating the lookup table took about ", elapsed_time) #

# peri.open_lookup("results/hydrus_loam")


# #create the global lookup table
# #t = time.time()
# #peri.create_lookup_global("results/"+peri.water_filename, sp)
# #elapsed_time = time.time() - t
# #print("Creating the global lookup table took about ", elapsed_time) #

# peri.open_global_lookup("results/"+peri.water_filename)

# # create lookup tables for the solute flow
# Ds_=np.logspace(np.log10(1.0e0), np.log10(1.0e1), 2)
# #peri.create_integralDiffusion_lookup("results/hydrus_loam_ss_solutes", sp)
# #peri.create_integralconcentration_lookup("results/hydrus_loam_sr_solutes", Ds_, sp)

# # open lookup tables for the solute flow
# peri.open_lookup_solutes("results/hydrus_loam_ss_solutes")
# peri.open_lookup_sr_solutes("results/hydrus_loam_sr_solutes")


# #numerically solve each root segment, this is the reference solution 
# hsr = np.array([peri.soil_root_interface_(rx[i], sx[i], inner_kr[i], rho[i], peri.sp) for i in range(0, len(rx))])
# tests.loc[:,"hsr"] = hsr
# print("norm of the root soil matrix potentials: ", LA.norm(hsr))

# #the standard implementation
# hsr1= peri.soil_root_interface_potentials_table(rx, sx, inner_kr, rho)
# tests.loc[:,"hsr_lookup"] = hsr1
# print("Norm of the difference to basic lookup table:", LA.norm(hsr-hsr1))

# #the alternative global implementation
# hsr3 = peri.soil_root_interface_potentials_table_global(rx, sx, inner_kr, rho)
# tests.loc[:,"hsr_global"] = hsr3
# print("Norm of the difference to global lookup table:", LA.norm(hsr-hsr3))

# #waterflow = 2*3.14*np.multiply((hsr-rx),inner_kr)/r_root
# waterflow = 2*3.14*np.multiply((sx-rx),inner_kr)#/r_root
# Phi_root = np.array([vg.fast_mfp[sp](hsr[i]) for i in range(len(hsr))])
# Phi_soil = np.array([vg.fast_mfp[sp](sx[i]) for i in range(len(sx))])
# tests.loc[:,"waterflow"] = waterflow
# Ds = 1.0e1

# #base solutes
# tests.loc[:,"sol_c"] = c_sol
# tests.loc[:,"sol_U"] = np.array([Vmax[i]*c_sol[i]/(Km[i]+c_sol[i])*1e6 for i in range(len(c_sol))])

# # steady state solutes
# solutes_ss = peri.soil_root_solutes_ss_(Phi_root, Phi_soil, c_sol, Vmax, Km, Ds, waterflow, peri.sp)
# tests.loc[:,"sol_c_ss"] = solutes_ss
# tests.loc[:,"sol_U_ss"] = np.array([Vmax[i]*solutes_ss[i]/(Km[i]+solutes_ss[i])*1e6 for i in range(len(solutes_ss))])

# # steady rate solutes
# solutes_sr = peri.soil_root_solutes_sr_(Phi_root, Phi_soil, rho, c_sol, Vmax, Km, Ds, waterflow, peri.sp)
# tests.loc[:,"sol_c_sr"] = solutes_sr
# tests.loc[:,"sol_U_sr"] = np.array([Vmax[i]*solutes_sr[i]/(Km[i]+solutes_sr[i])*1e6 for i in range(len(solutes_sr))])




# print("All inputs:")
# print(tests[inputs])
# print("All outputs:")
# print(tests[outputs])




