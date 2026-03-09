""" 
Q:
- Y less buds than leaves
- if we make kx age dependent better
- Y tested auxin and suc still at defualt val
"""

import types
import importlib
import os
import sys
import numpy as np
import psutil
import time
from pathlib import Path

SRC_PATH = "../../src/"
sys.path.append("../.."); sys.path.append(SRC_PATH)

# Create a fake plantbox namespace
plantbox = types.SimpleNamespace()

# Automatically import all folders inside src and attach to plantbox
for name in os.listdir(SRC_PATH):
    folder_path = os.path.join(SRC_PATH, name)
    if os.path.isdir(folder_path) and not name.startswith('__'):
        try:
            module = importlib.import_module(name)
            setattr(plantbox, name, module)
            sys.modules[f'plantbox.{name}'] = module
        except ModuleNotFoundError:
            # skip folders that are not importable as modules
            pass
            
import plantbox as pb  # |\label{l13:cplantbox}|
import plantbox.visualisation.vtk_plot as vp  # |\label{l13:vtk_plot}|
from plantbox.functional.phloem_flux import PhloemFluxPython 
from plantbox.functional.PlantHydraulicParameters import PlantHydraulicParameters 
sys.path.append("./modules"); 
from helpful import *

isCluster = (os.environ['HOME'] == '/home/m.giraud')


def runSim(directoryN_,nodeD, thread,   Qmax_, 
            thresholdSuc , BerthLim , doCondition, growthUpToNode = 9,simMax_ =  25 + 10/24 ): 
    
    # simulation type
    doDecapitation = (nodeD > 0)
    
    # for printing
    directoryN = directoryN_
    namef = 'results'+ directoryN+ "running"+ "_"+str(thread) + '.csv'
    my_file = Path(namef)
        
    # init printing
    write_file_float("running", thread, directoryN=directoryN)  # to know which thread are currently running
    
    
    # initialize vals    
    start_time = time.time()    
    outcondition = 0 # did the calibration fail? 
    temp_time  = time.time()
    weatherInit = weather(0, Qmax_)    
    simDuration = 0
    burninDuration = simDuration
    
    changedSimMax = False
    print('go to 1s for troubleshooting')
    dt = 1/24#60/60 #1h
    if simMax_ < 0:
        simMax = 100
    else:
        simMax = simMax_

    # plant system 
    pl = pb.MappedPlant(seednum = 2) 
    pl.readParameters("./modelparameter/UQ_1LeafRS.xml")
    depth = 60
    sdf = pb.SDF_PlantBox(np.inf, np.inf, depth )
    pl.setGeometry(sdf) # creates soil space to stop roots from growing out of the soil
    pl.initialize(verbose = False)
    
    ## todo: remove and set sucrose in seed instead?
    # makes plant old enough that it has enough leaf to grow on its own
    pl.simulate(3, False)
    simDuration += 3
        
        

    """ Coupling to soil """
    layers = depth; soilvolume = (depth / layers) * 3 * 12
    k_soil = []
    p_mean = weatherInit["p_mean"]#-187
    p_bot = p_mean + depth/2
    p_top = p_mean - depth/2
    sx = np.linspace(p_top, p_bot, depth)
    picker = lambda x,y,z : max(int(np.floor(-z)),-1) #abovegroud nodes get index -1
    pl.setSoilGrid(picker)  # maps segment

    """ Parameters phloem and photosynthesis """
    # TODO: add teh water stuff in there
    params = PlantHydraulicParameters()  # |\label{l52:hydraulic}|
    #setKrKx_xylem(20.,params)
    params.read_parameters("./modelparameter/pea_toredo")  # |\label{l52:hydraulic_end}|
    r = PhloemFluxPython(pl, params) #kss=0.2,kaa=1
    r.initialize()

    #r.wilting_point = -10000
    r.cpb_2_pm.maxLBud =  np.array([1.,1.,1.,1.])
    r.cpb_2_pm.maxLBudDormant = np.array([0.05,0.05,0.15,0.05])
    r.cpb_2_pm.budGR = 1.8 #same as fast growth for now #0.1
    r.cpb_2_pm.CSTthreshold = thresholdSuc  # sucrose threshold to release bud
    r.cpb_2_pm.L_dead_threshold = 2000. # if berthlim below this, bud dies
    r.cpb_2_pm.BerthLim = BerthLim # limit to go from bud to branch
    r.cpb_2_pm.limLenActive =  0.85 # limite size of 'young expanding' leaves
    r.auxin_D = 0. #muM /d   # auxin deletion rate
    r.auxin_P =  6.8e-3 #muM /d  # auxin production rate
    r.auxin_alpha =0.00239336
    r.initValAuxin = 0
    r.burnInTime = True
    
    #setKrKx_xylem(20., 0., r) 
    setKrKx_phloem(r.cpb_2_pm)
    setProtosynthesis_data(r)
    setPhloemflow_data(r) 
    r.cs = [weatherInit["cs"]]
    
    print('Rmax_st',r.cpb_2_pm.Rmax_st)
    
    step_id = -1
    while  simDuration < simMax:#
        step_id += 1
        
        temp_time = time.time()
        
        ####
        #   printing
        ####
        numLNodes = r.plant.getOrgans(3, False)[0].getNumberOfLaterals()
        Mstem = r.plant.getOrgans(3, False)[0]
        kids4distbase = np.array([Mstem.getChild(nkdb) for nkdb in range(Mstem.getNumberOfChildren())])
        
        MstemKidOT = np.array([kkk.organType() for kkk in kids4distbase])
        allLeaves = kids4distbase[MstemKidOT==pb.leaf]
        lastleafId = max(0,min(len(allLeaves),growthUpToNode)-1)
        lastLeaf = allLeaves[lastleafId]
        
        leafRank = len(allLeaves)
        
        leafArea =lastLeaf.getLength(True) * lastLeaf.getParameter("Width_blade")  
        maxLeafArea = lastLeaf.getParameter("k") * lastLeaf.getParameter("Width_blade") 
        
        expandedLeaf = leafArea > 0.1*maxLeafArea
        if __name__ == '__main__':
            print("leaf areas",np.array([ll.getLength(True) * ll.getParameter("Width_blade") for ll in allLeaves]))
            print("last leaf",lastleafId,leafArea,maxLeafArea,'num leaves',leafRank)
            print("do decapitate?",numLNodes, growthUpToNode, (not r.burnInTime), nodeD, expandedLeaf , (not changedSimMax))
        forDecapitate = (numLNodes >= growthUpToNode) and (not r.burnInTime) and (nodeD !=0) and expandedLeaf and (not changedSimMax)
        
        
        if (not doDecapitation) and (not changedSimMax) and (nodeD ==0) and (numLNodes >=growthUpToNode) and expandedLeaf: 
            if simMax_ <0:
                simMax = simDuration + testTime
            changedSimMax = True
            dt = 30/(60*24) #1MIN
            
        if r.burnInTime:
            weatherX = weather(burninDuration, Qmax_)
        else:
            weatherX = weather(simDuration+burninDuration, Qmax_)
                       
      
        r.solve(sim_time = simDuration, 
                rsx=sx, 
                 ea= weatherX["ea"], es= weatherX["es"], PAR = weatherInit["Qlight"],
                 TairC = weatherX["TairC"])    
                 
        r.solve_phloem_flow(dt, simDuration = simDuration ,TairC =  ( weatherX["TairC"]  +273.15))
      
    
        Q_Auxin =  r.Q_Auxin
        AuxinSource = np.array(r.cpb_2_pm.AuxinSource)
        orgs = r.plant.getOrgans(3)
        mainstem = orgs[0]
        maintstemNodesId = np.array(mainstem.getNodeIds())
        Q_Auxin_stem  = Q_Auxin[maintstemNodesId]
        C_Auxin_stem      = r.C_Auxin_np[maintstemNodesId]
        
        orgs = orgs[1:]
        auxTested= np.array([org.auxTested for org in orgs])
        lengthth_org = np.array([org.getLength(False) for org in orgs])
        sucTested = np.array([org.sucTested for org in orgs])
        budStage = np.array([org.budStage for org in orgs]) 
        BerthFact = np.array([org.BerthFact for org in orgs]) 
        orgsLeaves = r.plant.getOrgans(4, True)
        AuxinSourceK = np.array([AuxinSource[org.getNodeId(0)] for org in orgsLeaves]) 
         
        
        if __name__ == '__main__':
            print("\n\n\n\t\tat ", int(np.floor(simDuration)),"d", int((simDuration%1)*24),"h",
                        np.round(np.mean(r.Qlight))*1e6,"mumol m-2 s-1",r.burnInTime,doDecapitation)
            r.printSucroseData()
            
            print("sucTested",sucTested)
            print("sucTested*",sucTested*lengthth_org)
            print("auxTested",auxTested)
            print("lengthth_org",lengthth_org)
            print("BerthFact",BerthFact)
            print("budStage",budStage)
            print("AuxinSource",AuxinSourceK)
            print("C_Auxin_stem",C_Auxin_stem)
            print("sum aux source",sum(AuxinSourceK), sum(AuxinSource))
        
            
        if(changedSimMax):
            timeSinceDecap = simDuration - (simMax - testTime)
        else:
            timeSinceDecap = -1
        
        real_time = time.time()
        outcondition = doCondition(r,timeSinceDecap, thread,(real_time - start_time)/(60*60*24), 
                                    outcondition, nodeD)
        if (((temp_time - start_time)/(60*60*24) > 2) or (outcondition != 0) or (max(r.C_ST_np)>3)): #success or falur
            print('failure', temp_time - start_time, outcondition, max(r.C_ST_np))
            return -1
            
        
        if not r.burnInTime:
            r.plant.simulate(dt, False)
            simDuration += dt
        else:
            if burninDuration == 0.:                    
                mainStemAuxBU = np.array(C_Auxin_stem[1:])
                mainStemAuxBU[:] = -1
                
            burninDuration += dt
            mainStemAux_std = np.max(np.abs((C_Auxin_stem[1:]-mainStemAuxBU)))
            mainStemAuxBU = C_Auxin_stem[1:]
            mainStemAux_mean = np.mean(C_Auxin_stem[1:])
            if __name__ == '__main__':
                print("mainStemAux_std",mainStemAux_std, mainStemAux_mean)
            if mainStemAux_std < mainStemAux_mean*0.1:
                r.burnInTime = False
                
                raise Exception
                
    delete_file("running")
        


    
if __name__ == '__main__':
    start_time_ = time.time()   
    arg1 = "WT";arg2 = "high"
    nameSim = arg1
    lightLevel = arg2
    main_dir=os.environ['PWD']#dir of the file
    directoryN = "/"+os.path.basename(__file__)[:-3]+"_"+str(nameSim)+"_"+str(lightLevel)+"/"
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
    
    assert nameSim in np.array(["DW","SLM","WT","WTD"])
    assert lightLevel in np.array(["low","medium","high"])
    Qinit = 500*1e-6 
    GrRatio_ = 10
    GrRatioLats_ = 5
    GrRatioLeaf = 10
    if lightLevel == "medium":
        Qinit *= 0.86 
    elif lightLevel == "low":
        Qinit *= 0.64
        
    if nameSim == "DW":
        fileparam_ = "UQ_1Leaf_DW"
        Qmax__ = Qinit*0.5 #- 25e-6 #175*1e-6
        GrRatio_ = GrRatio_/3
        GrRatioLats_ = GrRatioLats_/3
        GrRatioLeaf = GrRatioLeaf/3
    else:
        fileparam_ = "UQ_1LeafRS"
        Qmax__ = Qinit
        
    if nameSim == "SLM":
        BerthLim_ = 1000
        thresholdSuc_ = 0
    else:
        BerthLim_ =  3.
        thresholdSuc_ = 3e-2 #*2*2
        
    if nameSim == "WTD":
        nodeD_ = 7
        simMax__ = -1
    else:
        nodeD_ = 0
        simMax__ = 25 + 10/24
        
    runSim(directoryN_ = directoryN, Qmax_ = Qmax__, 
            nodeD = nodeD_, thread = 0, thresholdSuc = thresholdSuc_, BerthLim = BerthLim_,
            doCondition = doConditionDefault)
            
    end_time_ = time.time()