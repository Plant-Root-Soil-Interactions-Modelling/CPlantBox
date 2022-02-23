"""plant example

adapted from "example_plant_anim", "example10_nodalGrowth"
"CPlantBox_tutorial_FSPM2020.ipynb", "test_leafparameter"
We present here parameters specific for leaves and stems

"""
import sys ;sys.path.append("../../.."); sys.path.append("../../../src/python_modules")
import numpy as np;

import plantbox as pb
import vtk_plot as vp
import matplotlib.pyplot as plt

##parameters for example:
adaptSeed = False
adaptLeaf = False
adaptStem = False
leafRadial = False #radial or not
anim = False
zoomLeafShape = True
export = False
getImage = True

##create plant:
plant = pb.Plant()
# Open plant and root parameter from a file
path = "../../../modelparameter/plantWithLeaves/" 
name ="Crypsis_aculeata_Clausnitzer_1994"
#other options:
#"Zea_mays_4_Leitner_2014"
#"Lupinus_albus_Leitner_2014"
#"Heliantus_Pagès_2013"
#"Crypsis_aculeata_Clausnitzer_1994"
#"Brassica_napus_a_Leitner_2010"
#"0"
#"example1e"
plant.readParameters(path + name + ".xml")

        

    
    
plant.initialize()

if anim:
    dt = 0.1
    N_ = 200
    min_ = np.array([0, -20, 0])/2
    max_ = np.array([20, 20, 30.])/2
    anim = vp.AnimateRoots(plant)
    anim.min = min_
    anim.max = max_
    anim.res = [1, 1, 1]
    anim.file = "results/example_plant"
    anim.avi_name = "results/example_"
    anim.plant = True
    anim.start()
    for i in range(0, N_):
        plant.simulate(dt, False)
        anim.root_name = "subType"
        anim.update()

if getImage:
    # Simulate
    if not anim:
        plant.simulate(30, True)
    # Plot, using vtk
    vp.plot_plant(plant, "organType")
    # zoom on leaf--theory--2D
    print("2D leaf shape of a full grown leaf")
    lorg = plant.getOrgans(pb.leaf)[1]
    lrp = lorg.getLeafRandomParameter()    
    leafRadial = (lrp.parametrisationType == 0)
    if leafRadial:
        N = len(lrp.leafGeometry)
        yy = np.linspace(0, lorg.leafLength(), N)
        geom_x, geom_y = [],[]
        for i, x in enumerate(lrp.leafGeometry):
            geom_x.extend(x)
            geom_y.extend([yy[i]] * len(x))
        geom_x = np.array(geom_x)
        geom_y = np.array(geom_y)        
        a  = lorg.leafArea() / lorg.leafLength() # scale radius
        plt.plot(geom_x * a, geom_y, "g*")
        plt.plot(-geom_x * a, geom_y, "g*")

    else:
        geom_x_a =  np.array([0])
        geom_x_b = np.array([ x[-1] for x in lrp.leafGeometry]) #normalized x value along length
        geom_x = np.concatenate((geom_x_a,geom_x_b))
        geom_y_a = np.array([0])
        geom_y_b =np.linspace(lrp.lb, lorg.leafLength()+lrp.lb, len(geom_x_b))
        geom_y = np.concatenate((geom_y_a,geom_y_b))
        a  = lorg.leafArea() / lorg.leafLength() # scale radius
        plt.plot(geom_x * a, geom_y, "g-*")
        plt.plot(-geom_x * a, geom_y, "g-*")
    plt.ylim([0, lrp.lmax+1])
    plt.xlim([-a-1, a+1])
    plt.axis('scaled')
    plt.show()
    
    # zoom on leaf--realized
    print("3D leaf shape of actual leaf")
    vp.plot_leaf(lorg)
    
    
