import os
#sourcedir = os.getcwd()+"/../../../"
sourcedir = "../../../CPlantBox/"
import sys;  
sys.path.append(sourcedir+"../dumux-rosi/python/modules");
sys.path.append(sourcedir+"../dumux-rosi/build-cmake/cpp/python_binding/");
sys.path.append(sourcedir)  
sys.path.append(sourcedir+"src")  
import matplotlib.pyplot as plt  
import numpy as np  
 
import plantbox as pb
import visualisation.vtk_plot as vp # for quick vizualisations
from rosi_richards import RichardsSPnum  # C++ part (Dumux binding)
from rosi_richards import RichardsSP
from richards import RichardsWrapper  # Python part  
from rosi_richards_cyl import RichardsCylFoam  # C++ part (Dumux binding), macroscopic soil model
from richards_flat import RichardsFlatWrapper  # Python part of cylindrical
import functional.van_genuchten as vg

def visualize_soil(soil, layers_ID, layers_pos):
    min_b = [-1., -1., -150.] # size [cm]  
    max_b = [0., 0., 0.] # size [cm]  
    cell_number = [1, 1, 150] # resolution
    s = RichardsWrapper(RichardsSPnum()) 
    s.initialize()     
    s.createGrid( min_b, max_b, cell_number)      
    s.setHomogeneousIC(-400., equilibrium = True)  # cm pressure head   
    s.setVGParameters(soil)
    s.setTopBC("noFlux")  
    s.setBotBC("noFlux")
    # define soil layers
    s.setLayersZ(layers_ID, layers_pos)

    s.initializeProblem()      
    
    s.solve(1/24/3600)
    points = s.getDofCoordinates()
    theta = s.getWaterContent() 
    print('ok')

# Define Mualem van Genuchten parameters for Selhausen soil profile according to Bauer et al. (2011, table 3, \url{https://doi.org/10.1007/s10533-011-9583-1}) |\label{l61ies:genuchten_a}|
# theta_r (-), theta_s (-), alpha (1/cm), n (-), Ks (cm d-1)
sand =  [0.045, 0.43, 0.15, 3, 1000]

l1 = [0.008, 0.389, 0.012, 1.97, 91.68]  # 0-20 cm 
l2 = [0.008, 0.389, 0.023, 1.23, 63.36]  # 20-33 cm
l3 = [0.008, 0.389, 0.01, 1.1, 10]  # 33-57 cm
l4 = [0.008, 0.389, 0.01, 1.1, 10]  # 57-120 cm
soil = [l1, l2, l3, l4] # Combine the hydraulic conductivity vectors from all soil layers to define soil type for simulation  

layers_ID = [4, 4, 3, 3, 2, 2, 1, 1]  
layers_pos = [-120., -57., -57., -33., -33, -20, -20, 0] 

visualize_soil(soil, layers_ID, layers_pos)