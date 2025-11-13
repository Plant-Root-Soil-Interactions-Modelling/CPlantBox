"""
Get the new parameters 
"""
import sys; 
from pathlib import Path
home = str(Path.home())
#edit CPlantBox dir
CPlantBox_dir = home + "/CPBLukas/CPlantBox" # "/data2model_0807/CPlantBox"

sys.path.append(CPlantBox_dir)
sys.path.append( CPlantBox_dir+"/src")

# sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import visualisation.vtk_plot as vp

import numpy as np
from structural.Plant import PlantPython
import matplotlib.pyplot as plt
import copy
import os
import pickle
from scipy.stats import gaussian_kde
import time

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
    
    
def run_benchmark(genotype = 'B'): #llambdao, kko,
    start_time = time.time()
    output = []
    for rep in range(1):
        N = 50
        outputs_12 = {
                'num':[0.,0.,0.,0.,0.],
                'length':[0.,0.,0.,0.],
                'ratio':0
                }
        outputs_349 = {
                'num':[0.,0.,0.,0.,0.],
                'length':[0.,0.,0.,0.],
                'ratio':0
                }
        outputs_50 = {
                'num':[0.,0.,0.,0.,0.],
                'length':[0.,0.,0.,0.],
                'ratio':0,
                'kde_main':[0. for i in range(50)],
                'kde_sub':[0. for i in range(50)],
                'kde_subsub':[0. for i in range(50)]
                }
                
        outpouts_mean = {
            'year'+str(i + 1): copy.deepcopy(outputs_12) if i < 2 else copy.deepcopy(outputs_349) if i < 49 else copy.deepcopy(outputs_50)
            for i in range(N)
        }
        
        long_root_types = np.array([1,2,3,4,5])
        fine_root_types = np.array([6,7,8,9])
        subtypes = max(long_root_types)
        
        orders = {'main' : [2], 'sub' : [3], 'subsub' : [4,5]}
        
        soilSpace = pb.SDF_PlantContainer(np.inf, np.inf,  np.inf, True)  # to avoid root growing aboveground
        plant = PlantPython()

        # Open plant and root parameter from a file
        plant.readParameters( CPlantBox_dir + '/experimental/wine_fichtl/' + genotype + "-wineV2.xml")

        plant.setGeometry(soilSpace) 


        yr_to_BEDD = 1225
        k_No = 1
        for ii, pp in enumerate(plant.getOrganRandomParameter(2)):
            
            if ii == 0:      
                pp.ldelay  = 0*yr_to_BEDD
                pp.ldelays = yr_to_BEDD * 50.#params['ldelays0'] #200#*5
                pp.successorP = [[1.]]# [[params['successorP0']]]
                pp.successorNo = [10]#[params['successorNo0']] 
                pp.successorST = [[2]]
            elif ii == 1:      
                pp.ldelay  = 0*yr_to_BEDD
                pp.ldelays = yr_to_BEDD *0.5#yr_to_BEDD * params['ldelays0'] #200#*5
                pp.successorP = [[1]]#[[params['successorP0']]]
                pp.successorNo = [1] #[params['successorNo0']] 
                pp.successorST = [[2]]
                pp.k_survive = 1.261
                pp.lambda_survive = 9.426
                pp.rlt_winter_max = 26.8 
                pp.rlt_winter_min =  3.169
            elif ii == 2:      
                pp.ldelay  = 0*yr_to_BEDD
                pp.ldelays = 0.5*yr_to_BEDD #yr_to_BEDD * params['ldelays0'] #200#*5
                pp.successorNo = [1 * k_No] #[params['successorNo0']] 
                pp.successorP = [[1./pp.successorNo[0], 0.], [0.1/pp.successorNo[0], 0.9]]#[[params['successorP0']]]
                pp.successorST = [[ii + 1, ii + 4], [ii + 1, ii + 4]]
                pp.successorP_age = [1225*3, 1225 * 100]
                pp.k_survive = 1.261
                pp.lambda_survive = 9.426
                pp.rlt_winter_max = 26.8 
                pp.rlt_winter_min =  3.169
            elif ii == 3:      
                pp.ldelay  = 0.#yr_to_BEDD
                pp.ldelays = yr_to_BEDD* 0.5#yr_to_BEDD * params['ldelays0'] #200#*5
                pp.successorNo = [1 * k_No] #[params['successorNo0']] 
                pp.successorP = [[0.7/pp.successorNo[0], 0.3/pp.successorNo[0]], [0.1/pp.successorNo[0], 0.9]]#[[params['successorP0']]]
                pp.successorST = [[ii + 1, ii + 4], [ii + 1, ii + 4]]
                pp.successorP_age = [1225*3, 1225 * 100]
                # pp.r *= 2
                pp.k_survive = 1.261
                pp.lambda_survive = 9.426
                #pp.k_survive = 1.367 
                #pp.lambda_survive =  4.995 
                pp.rlt_winter_max = 13
                pp.rlt_winter_min =  0
            elif ii == 4:      
                pp.ldelay  = 0*yr_to_BEDD
                pp.ldelays = yr_to_BEDD*0.5 #yr_to_BEDD * params['ldelays0'] #200#*5
                pp.successorNo = [1 * k_No] #[params['successorNo0']] 
                pp.successorP = [[0.5/pp.successorNo[0], 0.5/pp.successorNo[0]], [0.1/pp.successorNo[0], 0.9]]#[[params['successorP0']]]
                pp.successorST = [[ii + 1, ii + 4], [ii + 1, ii + 4]]
                #pp.k_survive = 1.261
                #pp.lambda_survive = 1.367
                #pp.rlt_winter_max = 13.9 
                #pp.rlt_winter_min =  0.
                pp.k_survive = 1.261
                pp.lambda_survive = 9.426
                pp.rlt_winter_max = 13/2
                pp.rlt_winter_min =  0
            elif ii == 5:      
                pp.ldelay  = 0*yr_to_BEDD
                pp.ldelays =yr_to_BEDD*0.5 #yr_to_BEDD * params['ldelays0'] #200#*5
                pp.successorP = [[1]]#[[params['successorP0']]]
                pp.successorNo = [1 * k_No] #[params['successorNo0']] 
                pp.successorST = [[ii + 4]]
                #pp.k_survive = 1.261
                #pp.lambda_survive = 1.367
                #pp.rlt_winter_max = 13.9 
                #pp.rlt_winter_min =  0.
                pp.k_survive = 1.261
                pp.lambda_survive = 9.426
                pp.rlt_winter_max = 13/2/2
                pp.rlt_winter_min =  0
            else: # fine roots
                pp.r = 1. # don t know the value, i just want something high
pp.lmax =                 
            
            
        '''
        start simulation
        '''
        # select randomly one of the rsml files of year 1
        file_path =  CPlantBox_dir + '/experimental/wine_fichtl/rsml/RSML_year1/'
        file_names = [name for name in os.listdir(file_path) if name[0] == genotype ] # repetitions of same genotype

        data_file = file_names[rep]
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
        # all_alive_long = []
        # all_alive_short = []
        dt = yr_to_BEDD  # ~1 yr

        # get3Dshape(plant,title_ = 'wine0', saveOnly = True)
        for i in range(N):
            print('age', i, end=", ", flush = True)
            plant.survivalTest()
            plant.simulate(dt, False) # i * 
            
            orgs_ = []
            for oo in plant.getOrgans(2, False):
                distance = 0.
                ooc = oo
                stc = ooc.param().subType
                while stc > 2:
                    distance += ooc.getParent().getLength(ooc.parentNI)
                    ooc = ooc.getParent()
                    stc = ooc.param().subType
                if distance < 175:
                    orgs_.append(oo)      
            
            all_ages = np.array([org.getAge()/yr_to_BEDD for org in orgs_ if org.isAlive()])
            all_alive.append(np.array([ org.isAlive() for org in orgs_  if  org.isAlive()]))
            all_subtypes.append(np.array([org.param().subType for org in orgs_  if  org.isAlive()]))
            all_lengths.append(np.array([org.getLength() for org in orgs_  if  org.isAlive()]))
            '''
            Ratio
            ''' 
            len_fine = 0
            is_fine_roots = np.isin(all_subtypes[i],fine_root_types)
            if sum(is_fine_roots) > 0:
                len_fine = sum(all_lengths[i][is_fine_roots])
            len_long = sum(all_lengths[i][np.invert(is_fine_roots)]) # ignore the length of 1 as seem different between xml and xlsx
            value = len_fine/(len_long + len_fine)
            print('percent', np.round(value*100), end=", ", flush = True)    
            '''
            Num
            '''
            nums = []
            for st in range(subtypes): # from 2 to 5            
                nums.append(sum(all_alive[i][all_subtypes[i] == (st + 1)]))
            nums.append(sum(all_alive[i][all_subtypes[i] > (st + 1)]))
            print('num', nums, end=", ", flush = True)
            # all_alive_long.append(np.array([org.isAlive() for org in orgs_ if org.param().subType in long_root_types]))
            # all_alive_short.append(np.array([org.isAlive() for org in orgs_ if org.param().subType in fine_root_types]))

            # # get3Dshape(plant,title_ = 'wine'+str(i+1), saveOnly = True) 
            
        # print('postprocessing')
        print("--- %s seconds for plant development---" % (time.time() - start_time),rep)
        
        for year in range(N):
            '''
            Lengths
            '''
            for st in range(1,subtypes): # from 2 to 5
                value = sum(all_lengths[year][all_subtypes[year] == (st + 1)])
                outpouts_mean['year'+str(year+1)]['length'][st-1] = value  
            '''
            Num
            '''
            for st in range(subtypes): # from 2 to 5            
                value = sum(all_alive[year][all_subtypes[year] == (st + 1)])
                outpouts_mean['year'+str(year+1)]['num'][st] = value
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
        idy = N - 1
        
        with open(CPlantBox_dir + '/experimental/wine_fichtl/measurements'+ genotype +'InitXX.pkl','rb') as f:
            xx = pickle.load(f)
             
        for key, value in orders.items():
            # List to store KDE y-values for each FileName
            
            x = xx[key]
            
            root_age = all_ages[np.isin(all_subtypes[idy],value)]
            
            if len(root_age) > 0:
                kde = gaussian_kde(root_age)
                y = kde(x)
            elif len(root_age) == 1:
                root_age.append(root_age[0])
            else:
                y = np.zeros(len(x))
            
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
        output.append(outpouts_mean)#list(flatten_values(outpouts_mean)))
    
        
        print("--- %s seconds for plant development---" % (time.time() - start_time),rep)
        
    #if getdict:
    #    return outpouts_mean
        
    # print('outpouts_mean',outpouts_mean)
    #f = [i for s in output for i in s]
    return output #np.array(f)
        


#return the filled parameter file from the array of parameter to test
def getFilledFile(X):
    
    params = {
        'ldelays0':X[0],'successorNo0':X[1].astype('int'),'successorP0' : X[2],
        'ldelays1':X[3],'successorNo1':X[4].astype('int') ,'successorP10':X[5]}
        
    return params

#return the filled parameter file from the array of parameter to test
def getDefaultFile():
    
    params = {
        'ldelays0':50,'successorNo0':10,'successorP0' : 1.,
        'ldelays1':40,'successorNo1':1 ,'successorP10':0.1
        }
        
    return params
    
if __name__ == '__main__': 
    # orders =  ['main', 'sub', 'subsub']
    # outpout_kde_mean = {}
    # outpout_kde_sd = {}
    # xx = {}
    # outpout_mean = {}
    # for genotype in ["B", "D", "E"]:

        # with open(CPlantBox_dir + '/experimental/wine_fichtl/measurements'+ genotype +'Init.pkl','rb') as f:
            # temp = pickle.load(f)
            # outpout_mean[genotype] = temp['mean']
    # for genotype in ["B", "D", "E"]:
        # outpout_kde_mean[genotype] = {}
        # for oo in orders:        
            # outpout_kde_mean[genotype][oo] = outpout_mean[genotype]['year50']['kde_'+oo]
            
        # with open( CPlantBox_dir + '/experimental/wine_fichtl/measurements'+ genotype +'InitXX.pkl','rb') as f:
                # temp = pickle.load(f)
                # xx[genotype] = temp
    # #params = {
    # #    'ldelays0':50,'successorNo0':10,'successorP0' : 1.,
    # #    'ldelays1':40,'successorNo1':1 ,'successorP10':0.2} #might be best to set No to 1 and change the ln
    # def getkde( lambdao, ko)
        # output = run_benchmark(getDefaultFile(), genotype = 'B', getdict = True, llambdao = lambdao, kko = ko)
        # return output['year50']['kde_sub']
    # #print(output)
    # mean_S = outpout_kde_mean['B']['kde_sub']
    # def target_length(x:np.array):
        # out = getkde(x[0], x[1])
        # sum = 0.
        # for i in range(len(mean_S)):
            # sum += (out[i] - mean_S[i]) ** 2
        # return np.sqrt(sum / (np.max([len(mean_S) - 2, 1])))  # -(k+1)
    # x0 =[4.995e+00,  1.367e+00]
            
    # f = lambda x0: target_length(x0) 
    # res = minimize(f, x0, method = 'Nelder-Mead', tol = 1e-6)
    output = run_benchmark(genotype = 'B')
    ratios = [output[0]['year'+str(year+1)]['ratio'] for year in range(50)]
    ages = [i+1 for i in range(50)]
    plt.plot(ages, ratios)
    plt.show()