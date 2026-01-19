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
from scipy.optimize import minimize

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
                
"""
pick path and files 
"""
def create_xml(genotype):
    file_names = []
    measurement_times = []
    yr_to_BEDD = 1225
    for year in [1,2]:
        file_path = "./rsml/RSML_year"+str(year)+"/"
        file_names_ = ["RSML_year"+str(year)+"/"+name for name in os.listdir(file_path) if name[0] == genotype ]
        file_names.extend(file_names_) # repetitions of same genotype
        measurement_times.extend([yr_to_BEDD*year for _ in range(len(file_names_))])  # days
        
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
        set fine roots
    """
    indices_type3b = [len(st3) for st3 in data.pick_order(3) ]
    indices_type4b = [len(st3) for st3 in data.pick_order(4) ]
    indices_type5b = [len(st3) for st3 in data.pick_order(5) ]
    
    for i, d in enumerate(data.rsmls):
        for j, _ in enumerate(d.polylines):
            length = data.rsmls[i].properties["length"][j]
            data.rsmls[i].properties["age"] = measurement_times[i]
            age = measurement_times[i]
            diam = np.mean(data.rsmls[i].functions["diameter"][j])
            assert ((age == yr_to_BEDD) or (age == yr_to_BEDD * 2))
            fry1 = (age == yr_to_BEDD) & (length <= 5.) & ( diam < 0.093)
            fry2 = (age == yr_to_BEDD * 2) & (length <= 5.) 
            if fry1 or fry2:
                data.rsmls[i].properties["order"][j] = data.rsmls[i].properties["order"][j] + 4
    indices_type3a = [len(st3) for st3 in data.pick_order(3) ]
    indices_type4a = [len(st3) for st3 in data.pick_order(4) ]
    indices_type5a = [len(st3) for st3 in data.pick_order(5) ]
    #print(indices_type3b)
    print(genotype)
    print('root st 3, year 1',np.mean(indices_type3a[:8]),'year 2',np.mean(indices_type3a[8:]))
    print('root st 4, year 1',np.mean(indices_type4a[:8]),'year 2',np.mean(indices_type4a[8:]))
    print('root st 5, year 1',np.mean(indices_type5a[:8]),'year 2',np.mean(indices_type5a[8:]))
    
            
    """ 
        find obvious parameters 
    """
    for i in range(0, max(subsub)):
        indices = data.pick_order(i)
        data.estimate_zones_(indices)  #  creates lb, ln, la, radius a, inseriton angle theta; e.g. writes single values into self.estimates[i][(j, "la")], where j is the root index
        data.aggregate_parameters_(indices, target_type = i)  # aggregates the individual root parameters (mean, sd) into data.parameters (list of RootRandomParameters) at index target_type

    """ 
        get growth rates         
    """
    st = 2
    indices_type2 = data.pick_order(st) 
    
    length_ = [[],[]]
    radii_ = [[],[]]
    for i, j_ in enumerate(indices_type2):
        for j in j_:
            year = int(data.rsmls[i].properties["age"]  > yr_to_BEDD)
            length_[year].append(data.rsmls[i].properties["length"][j])
            radii_[year].append(np.mean(data.rsmls[i].functions["diameter"][j]))
    radii_ = [np.array(radii_[0]),np.array(radii_[1])] 
    length_ = [np.array(length_[0]),np.array(length_[1])] 
    
    
    age_th = []
    for did, dd in enumerate(radii_[0])  :
        age_th.append((dd )/0.083 ) #- min(radii_[0])  * (365/170)
    age_th = np.array(age_th)
    
    tokeep = (age_th > 0) #& (age_th < yr_to_BEDD)
    age_scaled = (radii_[0]) / (radii_[0].max()) * yr_to_BEDD
    
    r_4_min = min(length_[0][age_scaled == min(age_scaled)]/ age_scaled[age_scaled == min(age_scaled)]) * 0.8

    # if st == 2:
        # def negexp_length(t, r, k):
            # return k * (1 - np.exp(-(r / max(k, 1.e-9)) * t))
          
        # # we set yr_to_BEDD * 1.5 for st 2 because almost all the st 2 seem born the 1st years (0.5) 
        # # then go through the 2nd yea (1)
        # plt.scatter(age_scaled,length_[0], label='obs_yr1')
        # plt.scatter([yr_to_BEDD * 1.5 for i in range(len(length_[1])) ],length_[1], label='obs_yr2')
        # plt.scatter(age_scaled, [negexp_length(t, r_4_min, lmax[genotype][2] ) for t in age_scaled], label='th')
        # #plt.scatter([yr_to_BEDD  for i in range(len(length_[1]))], length_[1], label = 'obs_yr2')
        # plt.scatter([yr_to_BEDD *  1.5 ], [negexp_length(yr_to_BEDD * 1.5, r_4_min,  lmax[genotype][2])], label = 'th_yr2')
        # plt.legend()
        # plt.show()

        # #raise Exception
    st = 3
    
    indices_type2 = data.pick_order(st) 
    
    length_ = [[],[]]
    for i, j_ in enumerate(indices_type2):
        for j in j_:
            year = int(data.rsmls[i].properties["age"]  > yr_to_BEDD)
            length_[year].append(data.rsmls[i].properties["length"][j])
    length_ = [np.array(length_[0]),np.array(length_[1])] 
    length_f = [item for sublist in length_ for item in sublist] 
    age_ = [[yr_to_BEDD*0.5 for i in range(len(length_[0]))], [yr_to_BEDD*1.5 for i in range(len(length_[1]))]]
    age_f = [item for sublist in age_ for item in sublist] 
    r_4_min_st3 = np.mean(length_[0])/(yr_to_BEDD*0.5)

    
    def negexp_length(t, r):
        return lmax[genotype][st] * (1 - np.exp(-(r / max(lmax[genotype][st], 1.e-9)) * t))
        
    def target_length(x:np.array):
        out = [negexp_length(yy, x[0]) for yy in age_f]
        sum = 0.
        for i in range(len(length_f)):
            sum += (out[i] - length_f[i]) ** 2
        return np.sqrt(sum / (np.max([len(age_f) - 2, 1])))  # -(k+1)

    f = lambda x0: target_length(x0) 
    x0 = [r_4_min_st3]
    # res = minimize(f, x0, method = 'Nelder-Mead', tol = 1e-6)
    
    # r_4_min_st3 = res.x[0]

    # if st == 3:
        # plt.scatter([yr_to_BEDD*0.5 for i in range(len(length_[0])) ],length_[0], label='obs_yr1')
        # plt.scatter([yr_to_BEDD*1.5 for i in range(len(length_[1])) ],length_[1], label='obs_yr2')
        # plt.plot([yr_to_BEDD*0.5, yr_to_BEDD*1.5],[negexp_length(t, r_4_min_st3) for t in [yr_to_BEDD*0.5, yr_to_BEDD*1.5]], label='th')
        # # plt.plot([yr_to_BEDD*0.5, yr_to_BEDD*1.5],[negexp_length(t, res.x[0]) for t in [yr_to_BEDD*0.5, yr_to_BEDD*1.5]], label='th2')
        # plt.legend()
        # plt.show()
    
    
        
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
    data.pparameters.Lmax_unsuberized = 4.
    data.pparameters.Lmax_suberized = 10.
    data.pparameters.delayDefinition = 4
    data.pparameters.seedPos = pb.Vector3d(0.,0.,-10.)

    
    for ii, pp in enumerate(data.parameters):
        pp.subType = ii
        pp.tropismS = 0.2
        pp.successorNo = [1]
        pp.a_s = 0.
        pp.la = max(pp.la, 0.)
        pp.dxMin = 0.2
        if ii <= max(subsub):
            pp.a = 0.093/2
            pp.a_gr =  0.083/2/yr_to_BEDD
            pp.is_fine_root = False
            if ii > 0 :
                pp.k_survive = survive[ii]['k_survive']
                pp.lambda_survive = survive[ii]['lambda_survive']
                pp.rlt_winter_max = survive[ii]['rlt_winter_max']
                
            if ii < main:
                pp.successorST =  np.array([[[main]]]) 
                pp.successorOT =  np.array([[[2]]]) 
                pp.successorP = np.array([[[1]]])
            else:
                pp.successorST =  np.array([ [[ii + 1]], [[ii + 4]] ]) 
                pp.successorOT =  np.array([ [[2]],       [[2]] ]) 
                pp.successorP = np.array([ [[1.]], [[0.]] ]) #0.64,0.35
                pp.lmax = lmax[genotype][ii]
                if ii == 2:
                    pp.r = r_4_min
                else:
                    pp.r = r_4_min_st3 # pp.lmax * growth_beta[ii]
                if ii == main:
                    pp.tropismT = 7
                    pp.tropismN = 0.5  
                    pp.tropismW1 = tropism[genotype][0] #0.85 #0.4 # gravitropism
                    pp.tropismW2 = tropism[genotype][1] # plagiotropism
            if ii == max(subsub):
                pp = data.parameters[min(subsub)].copy(plant)
                pp.subType = ii
                pp.successorST =  np.array([[[ii + 4]]]) 
                pp.successorOT =  np.array([[[2]]]) 
                pp.successorP = np.array([[[0.35]]])
        else:
            pp.a = 0.05/2
            pp.is_fine_root = True
            pp.lmax = 5.
            pp.r = data.parameters[max(subsub)].r
            pp.successorP = np.array([[[0.5]]])
            pp.successorST =  np.array([[[ii + 1]]]) 
            pp.successorOT =  np.array([[[2]]]) 
            # update?
            pp.la = data.parameters[min(subsub)].la
            pp.lb = data.parameters[min(subsub)].lb
            pp.ln = data.parameters[min(subsub)].ln
            pp.las = data.parameters[min(subsub)].las
            pp.lbs = data.parameters[min(subsub)].lbs
            pp.lns = data.parameters[min(subsub)].lns
            
        data.parameters[ii] = pp
            
    data.parameters[-1].successorST =  [] 
    
    """ write set """

    for ii, pp in enumerate(data.parameters):
        param = pp.copy(plant)
        plant.setOrganRandomParameter(param)
    
    param = data.pparameters.copy(plant)
    plant.setOrganRandomParameter(param)
    
    os.makedirs("results/xmlFiles", exist_ok=True)
    plant.writeParameters("./results/xmlFiles/"+genotype + "-wineV2.xml")
    return max(length_[0]),max(length_[1])
        
lb1, lb2 = create_xml('B')
le1, le2 = create_xml('E')
ld1, ld2 = create_xml('D')
print(le1, le2, ld1, ld2, lb1, lb2)