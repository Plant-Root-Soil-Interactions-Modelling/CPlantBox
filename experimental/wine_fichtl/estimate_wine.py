"""
specialized scipt for the wine rsml data 
"""
import sys; sys.path.append("../.."); sys.path.append("../../src/")
sys.path.append("../../gui/estimate/")
import visualisation.vtk_plot as vp
from estimate_data import EstimateDataModel
import estimate_plots as ep
import os
import matplotlib.pyplot as plt
import numpy as np
import plantbox as pb 

SMALL_SIZE = 16
MEDIUM_SIZE = 16
BIGGER_SIZE = 16
plt.rc('font', size = SMALL_SIZE)  # controls default text sizes
plt.rc('axes', titlesize = SMALL_SIZE)  # fontsize of the axes title
plt.rc('axes', labelsize = MEDIUM_SIZE)  # fontsize of the x and y labels
plt.rc('xtick', labelsize = SMALL_SIZE)  # fontsize of the tick labels
plt.rc('ytick', labelsize = SMALL_SIZE)  # fontsize of the tick labels
plt.rc('legend', fontsize = SMALL_SIZE)  # legend fontsize
plt.rc('figure', titlesize = BIGGER_SIZE)  # fontsize of the figure title
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

"""
pick path and files 
"""
def create_xml(genotype):
    file_names = []
    measurement_times = []
    for year in [1,2]:
        file_path = "./rsml/RSML_year"+str(year)+"/"
        file_names_ = ["RSML_year"+str(year)+"/"+name for name in os.listdir(file_path) if name[0] == genotype ]
        file_names.extend(file_names_) # repetitions of same genotype
        measurement_times.extend([365*year for _ in range(len(file_names_))])  # days
        
    """
    RSMLs:
    Type 0: Stem
    Type 1: the initially planted roots
    Type 2: first order roots => Main
    Type 3: second order roots => sub
    Type 4-5 => subsub
    Type 7-10: thin roots
    get: la, lb, theta, a (if available)
    """
    plant = pb.Plant()
    yr_to_BEDD = 1225
    numTypes = 10
    main = 2
    sub = 3
    subsub = [4, 5]
    """ 
        open rsml files 
    """
    data = EstimateDataModel(numTypes)  # new data model
    data.open_files("./rsml", file_names)  # open rmsl
    data.times = measurement_times  # initialize
    for i in range(0, len(data.times)):
        data.estimates[i] = {}

    """ 
        find obvious parameters 
    """
    for i in range(0, max(subsub)):
        indices = data.pick_order(i)
        data.estimate_zones_(indices)  #  creates lb, ln, la, radius a, inseriton angle theta; e.g. writes single values into self.estimates[i][(j, "la")], where j is the root index
        data.aggregate_parameters_(indices, target_type = i)  # aggregates the individual root parameters (mean, sd) into data.parameters (list of RootRandomParameters) at index target_type

        
        
    """
    from data:
    - max length
    - growth rate
    - diameter and diamter growth
    - k_survive, lambda_survive, tropismN, S, W1, W2, is_fine_root
    - successors
    + seed data
    """
    
    ## seed data
    data.pparameters.Lmax_unsuberized = 5.
    data.pparameters.Lmax_suberized = 10.
    data.pparameters.delayDefinition = 4

    ## root data
    growth_beta = {2: 3.28 * 10e-5,
                    3: 2.16 * 10e-5,
                    4: 2.00 * 10e-5,
                    5: 2.00 * 10e-5}
                    
    lmax = {'E':{}, 'B':{}, 'D':{}}
    lmax['E'][2] = 790
    lmax['E'][3] = 120
    lmax['E'][4] = 60
    lmax['E'][5] = 60

    lmax['B'][2] = 630
    lmax['B'][3] = 100
    lmax['B'][4] = 50
    lmax['B'][5] = 50

    lmax['D'][2] = 640
    lmax['D'][3] = 90
    lmax['D'][4] = 50
    lmax['D'][5] = 50
    
    tropism = {'E':[0.95,0.05], 'B':[0.85,0.15], 'D':[0.65,0.35]}
    
    survive = {1:{'k_survive':1.88,'lambda_survive':7.53,'rlt_winter_max':23},
                2:{'k_survive':1.88,'lambda_survive':7.53,'rlt_winter_max':23},
                3:{'k_survive':1.42,'lambda_survive':3.54,'rlt_winter_max':14}, 
                4:{'k_survive':2.17,'lambda_survive':4.55,'rlt_winter_max':8},
                5:{'k_survive':2.17,'lambda_survive':4.55,'rlt_winter_max':8} }
    
    for ii, pp in enumerate(data.parameters):
        pp.subType = ii
        pp.tropismS = 0.2
        pp.successorNo = [1]
        pp.a_s = 0.
        pp.la = max(pp.la, 0.)
        if ii <= max(subsub):
            pp.a = 0.093/2
            pp.a_gr =  0.083/2/yr_to_BEDD
            pp.is_fine_root = False
            if ii > 0 :
                pp.k_survive = survive[ii]['k_survive']
                pp.lambda_survive = survive[ii]['lambda_survive']
                pp.rlt_winter_max = survive[ii]['rlt_winter_max']
                
            if ii < main:
                pp.successorST =  np.array([[main]]) 
                pp.successorOT =  np.array([[2]]) 
                pp.successorP = np.array([[1]])
            else:
                pp.successorST =  np.array([[ii + 1, ii + 4]]) 
                pp.successorOT =  np.array([[2, 2]]) 
                pp.successorP = np.array([[0.64,0.35]])
                pp.lmax = lmax[genotype][ii]
                pp.r = pp.lmax * growth_beta[ii]
                if ii == main:
                    pp.tropismT = 7
                    pp.tropismN = 2  
                    pp.tropismW1 = tropism[genotype][0] #0.85 #0.4 # gravitropism
                    pp.tropismW2 = tropism[genotype][1] # plagiotropism
            if ii == max(subsub):
                pp = data.parameters[min(subsub)].copy(plant)
                pp.subType = ii
                pp.successorST =  np.array([[ii + 4]]) 
                pp.successorOT =  np.array([[2]]) 
                pp.successorP = np.array([[0.35]])
        else:
            pp.a = 0.05/2
            pp.is_fine_root = True
            pp.lmax = 5.
            pp.r = data.parameters[max(subsub)].r
            pp.successorP = np.array([[0.5]])
            pp.successorST =  np.array([[ii + 1]]) 
            pp.successorOT =  np.array([[2]]) 
            # update?
            pp.la = 1
            pp.lb = 1
            pp.ln = 1
        data.parameters[ii] = pp
            
    data.parameters[-1].successorST =  [] 
    
    """ write set """

    for ii, pp in enumerate(data.parameters):
        param = pp.copy(plant)
        plant.setOrganRandomParameter(param)
    
    param = data.pparameters.copy(plant)
    plant.setOrganRandomParameter(param)
    
    plant.writeParameters(genotype + "-wine.xml")
        

create_xml('B')
create_xml('D')
create_xml('E')