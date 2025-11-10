"""
Get the new parameters 
"""
import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import visualisation.vtk_plot as vp

import numpy as np
from structural.Plant import PlantPython
import matplotlib.pyplot as plt
import copy
import os
import pickle
from scipy.stats import gaussian_kde

def get3Dshape(plant,title_ = 'wine', saveOnly = True):        
    orgs_ = plant.getOrgans(2)
    orgs_a = [1 for org in orgs_ if org.isAlive() ]
    orgs_al = [1 for org in orgs_ if (org.isAlive() and not org.getParameter('is_fine_root')) ]
    ana = pb.SegmentAnalyser(plant) 
    segOs = plant.getSegmentOrigins()
    
    #vp.plot_roots(ana, "id")

    '''
    Lignification status, Survival, fine roots
    '''
    lignification = [segO.lignificationStatus() for segO in segOs]
    aliveSegs = [segO.isAlive() for segO in segOs]
    is_fine_root = [segO.getParameter('is_fine_root') for segO in segOs]
    
    print('alive nodes:', sum(aliveSegs), 'alive organs',sum(orgs_a), 'alive long lived roots',sum(orgs_al))
    # ana.addData('alive', aliveSegs)
    # ana.addData('lignification', lignification)
    # ana.addData('is_fine_root', is_fine_root)
    # ana.filter('alive', 1)
    # vp.plot_roots(ana, "subType",p_names = ['lignification','is_fine_root',"creationTime","id"] , win_title = title_, render = not saveOnly)
    
    
def run_benchmark(params, genotype = 'B'): 

    outputs_12 = {
            'num':[0.,0.,0.,0.,0.],
            'length':[0.,0.,0.,0.],
            'ratio':0
            }
    outputs_3to10 = {
            'num':[0.,0.,0.,0.,0.],
            'ratio':0
            }
    outputs_50 = {
            'num':[0.,0.,0.,0.,0.],
            'ratio':0,
            'kde_main':[0. for i in range(50)],
            'kde_sub':[0. for i in range(50)],
            'kde_subsub':[0. for i in range(50)]
            }
            
    outpouts_mean = {
        'year'+str(i+1): copy.deepcopy(outputs_12) if i < 2 else copy.deepcopy(outputs_3to10) if i < 49 else copy.deepcopy(outputs_50)
        for i in range(50)
    }
    
    long_root_types = np.array([1,2,3,4,5])
    fine_root_types = np.array([6,7,8,9])
    subtypes = max(long_root_types)
    
    orders = {'main' : [2], 'sub' : [3], 'subsub' : [4,5]}
    
    soilSpace = pb.SDF_PlantContainer(np.inf, np.inf,  np.inf, True)  # to avoid root growing aboveground
    plant = PlantPython(1)

    # Open plant and root parameter from a file
    plant.readParameters( genotype + "-wine.xml")

    plant.setGeometry(soilSpace) 


    yr_to_BEDD = 1225

    params['successorP11'] = 1. - params['successorP10']
    for ii, pp in enumerate(plant.getOrganRandomParameter(2)):
        
        if ii == 0:      
            pp.ldelay  = 0*yr_to_BEDD
            pp.ldelays = yr_to_BEDD * params['ldelays0'] #200#*5
            pp.successorP = [[params['successorP0']]]
            pp.successorNo = [params['successorNo0']] 
            pp.successorST = [[2]]
        else:
            pp.ldelay  = 0*yr_to_BEDD
            pp.ldelays = yr_to_BEDD * params['ldelays1']
            pp.successorNo = [params['successorNo1']] 
            if ii in np.array([5,6,7,8]):
                pp.successorP = [[params['successorP11']]]                 
            elif ii in np.array([2,3,4]):            
                pp.successorP = [[params['successorP10'], params['successorP11']]] 
            
        # print(ii,pp.successorST,pp.successorP,pp.lmax)
        
        
    '''
    start simulation
    '''
    # select randomly one of the rsml files of year 1
    file_path = "./rsml/RSML_year1/"
    file_names = [name for name in os.listdir(file_path) if name[0] == genotype ] # repetitions of same genotype

    data_file = file_names[int(np.round( plant.rand()*(len(file_names) - 1)))]
    plant.initialize_static(file_path + data_file, [0, 1])  # 0 is shoot, 1 are static roots

    # the static laterals 2 and 3 are replaced with the growing lateral 2
    ld, ld1 = plant.set_identical_laterals([0, 1], [1, 2, 3], 2)
    # plt.hist(np.array(ld1)/yr_to_BEDD, density = False, bins=30)
    # plt.title("Creation time of the main roots")
    # plt.show()
    # plt.hist(np.array(ld)/yr_to_BEDD, density = False, bins=30)
    # plt.title("Creation time of the main roots")
    # plt.show()

    plant.initialize_static_laterals()
    plant.betaN = 5000
      

    all_lengths = []
    all_ages = []
    all_subtypes = []
    all_alive = []
    all_alive_long = []
    all_alive_short = []
    dt = yr_to_BEDD  # ~1 yr
    N = 50

    # get3Dshape(plant,title_ = 'wine0', saveOnly = True)
    for i in range(N):
        # print('age',i)
        plant.survivalTest()
        plant.simulate(dt, False)
        orgs_ = plant.getOrgans(2)
        all_ages.append(np.array([org.getAge()/yr_to_BEDD for org in orgs_]))
        all_lengths.append(np.array([org.getLength() for org in orgs_]))
        all_alive_long.append(np.array([org.isAlive() for org in orgs_ if org.param().subType in long_root_types]))
        all_alive_short.append(np.array([org.isAlive() for org in orgs_ if org.param().subType in fine_root_types]))
        all_alive.append(np.array([org.isAlive() for org in orgs_]))
        all_subtypes.append(np.array([org.param().subType for org in orgs_]))
        # print(np.array([org.param().subType for org in orgs_]))
        # print(np.array([org.isAlive() for org in orgs_]))    
        # # get3Dshape(plant,title_ = 'wine'+str(i+1), saveOnly = True) 
    # # raise Exception
    # print('postprocessing')
            
    for year in range(2):            
        '''
        Lengths
        '''
        for st in range(1,subtypes): # from 2 to 5
            value = sum(all_lengths[year][all_subtypes[year] == (st + 1)])
            outpouts_mean['year'+str(year+1)]['length'][st-1] = value
            

    for year in range(N):
        '''
        Num
        '''
        for st in range(subtypes): # from 2 to 5
            value = sum(all_alive[year][all_subtypes[year] == (st + 1)])
            outpouts_mean['year'+str(year+1)]['num'][st] = value   
            # print(year, st,all_subtypes[year] == (st + 1),value )
            # #raise Exception
        # print(all_subtypes[year])
        # print(all_alive[year])
        # print('outpouts_mean["year"+str(year+1)]["num"]',outpouts_mean['year'+str(year+1)]['num'])
        # raise Exception
        '''
        Ratio
        ''' 
        len_fine = 0
        is_fine_roots = np.isin(all_subtypes[year],fine_root_types)
        if sum(is_fine_roots) > 0:
            len_fine = sum(all_lengths[year][is_fine_roots])
        len_long = sum(all_lengths[year][np.invert(is_fine_roots)]) # ignore the length of 1 as seem different between xml and xlsx
        value = len_fine/(len_long + len_fine)
        outpouts_mean['year'+str(year+1)]['ratio'] = value         
        

    '''
    Age distribution at last time step
    '''
    year = N - 1 
    
    with open('./measurements'+ genotype +'InitXX.pkl','rb') as f:
        xx = pickle.load(f)
         
    for key, value in orders.items():
        # List to store KDE y-values for each FileName
        
        x = xx[key]
        
        root_age = all_ages[year][np.isin(all_subtypes[year],value)]
        
        kde = gaussian_kde(root_age)
        y = kde(x)
        
        outpouts_mean['year50']['kde_' +key] = y

    def flatten_values(obj):
        if isinstance(obj, dict):
            for v in obj.values():
                yield from flatten_values(v)
        elif isinstance(obj, (list, tuple, np.ndarray)):
            for v in obj:
                yield from flatten_values(v)
        else:
            yield obj
    # print('outpouts_mean',outpouts_mean)
    return np.array(list(flatten_values(outpouts_mean)))   
        


#return the filled parameter file from the array of parameter to test
def getFilledFile(X):
    
    params = {
        'ldelays0':X[0],'successorNo0':X[1].astype('int'),'successorP0' : X[2],
        'ldelays1':X[3],'successorNo1':X[4].astype('int') ,'successorP10':X[5]}
        
    return params

#return the filled parameter file from the array of parameter to test
def getDefaultFile(X):
    
    params = {
        'ldelays0':50,'successorNo0':10,'successorP0' : 1.,
        'ldelays1':40,'successorNo1':10 ,'successorP10':0.5}
        
    return params
    
if __name__ == '__main__': 
            
    params = {
        'ldelays0':50,'successorNo0':10,'successorP0' : 1.,
        'ldelays1':40,'successorNo1':1 ,'successorP10':0.2} #might be best to set No to 1 and change the ln
    
    output = run_benchmark(params, genotype = 'B')
    print(len(output))