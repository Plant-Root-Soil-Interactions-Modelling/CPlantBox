""" water movement within the root (static soil) """
#directoryN = "/7to14dry/"
import sys; 
CPBdir = "../.."
sys.path.append(CPBdir+"/src");
sys.path.append(CPBdir);
sys.path.append("../../..");sys.path.append(".."); 
sys.path.append(CPBdir+"/src/python_modules");
sys.path.append("../build-cmake/cpp/python_binding/") # dumux python binding
sys.path.append("../../build-cmake/cpp/python_binding/")
sys.path.append("../modules/") # python wrappers 

#import matplotlib
#matplotlib.use('AGG') #otherwise "import matplotlib.pyplot" hangs

#from Leuning_speedup import Leuning #about 0.7 for both
#from photosynthesis_cython import Leuning
import plantbox as pb
#from plantbox import Photosynthesis as ph
#import vtk_plot as vp
import math
import os
import numpy as np
#import vtk_plot as vp
#import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import visualisation.vis_tools as cpbvis
import visualisation.vtk_plot as vp



isCluster = (os.environ['HOME'] == '/home/m.giraud')



    
""" Parameters """
def launchUQR(directoryN,simInit, condition,forPlants):
           
    
    simDuration = simInit # [day] init simtime
    simMax =simInit #+7
    depth = 60
    dt = 1/24 #1h
    verbose = True
    leaf_res = 30

    # plant system 
    pl = pb.MappedPlant(seednum = 2) #pb.MappedRootSystem() #pb.MappedPlant()
    pl2 = pb.MappedPlant(seednum = 3) #pb.MappedRootSystem() #pb.MappedPlant()
    path = CPBdir+"/modelparameter/structural/plant/"
    name = "Triticum_aestivum_test_2021"#"wheat_uqr15" #"manyleaves"## "Anagallis_femina_Leitner_2010"  # Zea_mays_1_Leitner_2010

    pl.readParameters(path + name + ".xml")
    pl2.readParameters(path + name + ".xml")
    
    if forPlants == "mix":
        for p in pl.getOrganRandomParameter(pb.root):
            #p.theta =  p.theta * 2
            #p.tropismS = p.tropismS/10
            p.tropismT = 0
    if forPlants == "deep":
        pass
    if forPlants == "shallow":
        for p in pl.getOrganRandomParameter(pb.root):
            p.theta =  p.theta * 2
            p.tropismS = p.tropismS/10
        for p in pl2.getOrganRandomParameter(pb.root):
            p.theta =  p.theta * 2
            p.tropismS = p.tropismS/10
                
    def setleaf(plant):            
        for p in plant.getOrganRandomParameter(pb.leaf):
          p.la,  p.lmax = 38.41053981, 38.41053981
          #p.theta = 0.2 # 0.2
          p.theta = 0.01
          p.tropismT = 1 # 6: Anti-gravitropism to gravitropism
          p.areaMax = 54.45388021  # cm2, area reached when length = lmax
          phi = np.array([-90,-80, -45, 0., 45, 90]) / 180. * np.pi    
          l = np.array([38.41053981,1 ,1, 0.3, 1, 38.41053981]) #distance from leaf center
          p.leafGeometryPhi = phi
          p.leafGeometryX = l
          #p.tropismN = 5
          p.tropismS = 0.08
          p.tropismAge = 5 #< age at which tropism switch occures, only used if p.tropismT = 6
          p.createLeafRadialGeometry(phi,l,leaf_res)
    
    setleaf(pl)
    setleaf(pl2)

    #raise Exception
    sdf = pb.SDF_PlantBox(np.inf, np.inf, depth )

    pl.setGeometry(sdf) # creates soil space to stop roots from growing out of the soil
    pl2.setGeometry(sdf) # creates soil space to stop roots from growing out of the soil


    pl.getOrganRandomParameter(pb.seed)[0].seedPos = pb.Vector3d(1.5,0, -0.3) 
    pl2.getOrganRandomParameter(pb.seed)[0].seedPos = pb.Vector3d(-1.5,0, -0.3) 
    
    pl.initialize(verbose = True, stochastic = True)
    pl2.initialize(verbose = True, stochastic = True)
    
    print("simulate")
    pl.simulate(simDuration, False)#, "outputpm15.txt")
    pl2.simulate(simDuration, False)#, "outputpm15.txt")
    
    """ Coupling to soil """
    

    """ Parameters phloem and photosynthesis """
    
    def doVis(plant, nameimage):    
        vis = pb.PlantVisualiser(plant)
        vis.SetGeometryResolution(8)
        vis.SetLeafResolution(leaf_res)
        vis.ResetGeometry()
        vis.ComputeGeometryForOrganType(pb.stem, False)
        vis.ComputeGeometryForOrganType(pb.leaf, False)
        vis.ComputeGeometryForOrganType(pb.root, False)

        # Write the geometry to file#
        data = cpbvis.PolydataFromPlantGeometry(vis)
        cpbvis.WritePolydataToFile(data,  'results' + directoryN+nameimage + ".vtp")

        #vp.plot_plant(plant, "subType")
        vp.write_plant( 'results' + directoryN+nameimage + "vtk",plant) # will write test.vtp (with 

    doVis(pl, "pl1")
    doVis(pl2, "pl2")
    print("simDuration", simDuration, "d")

    

if __name__ == '__main__':    
    directoryN = "/"+os.path.basename(__file__)[:-3]+"/"

    main_dir=os.environ['PWD']#dir of the file
    results_dir = main_dir +"/results"+directoryN
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)
    else:
        import shutil
        shutil.rmtree(results_dir)
        os.makedirs(results_dir)

    launchUQR(directoryN,7, "dry")