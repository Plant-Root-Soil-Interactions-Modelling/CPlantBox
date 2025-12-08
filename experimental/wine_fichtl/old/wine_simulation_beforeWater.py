"""
Run simulation
"""
import sys; 
from pathlib import Path

#edit CPlantBox dir
CPlantBox_dir =  "../.." # "/data2model_0807/CPlantBox"

from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank(); max_rank = comm.Get_size()
sys.path.append(CPlantBox_dir)
sys.path.append( CPlantBox_dir+"/src")
sys.path.append( CPlantBox_dir+"/gui/viewer")

# sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import visualisation.vtk_plot as vp
import viewer_conductivities

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
    
    ana.addData('lignification', lignification)
    ana.addData('is_fine_root', is_fine_root)
    ana.filter('alive', 1)
    vp.plot_roots(ana, "subType",p_names = ['lignification','is_fine_root',"creationTime","id"] , win_title = title_, render = not saveOnly)
    

def run_benchmark(xx, genotype = 'B', rep_input = -1): #llambdao, kko,
    start_time = time.time()
    output = []
    file_path =  CPlantBox_dir + '/experimental/wine_fichtl/rsml/RSML_year1/'
    file_names = [name for name in os.listdir(file_path) if name[0] == genotype ] # repetitions of same genotype

    data_file = file_names[rep_input]
    reps = 1# len(file_names)
    doVTP = False# (rep_input == 0)
    for rep in range(reps):
        N = 50
        outputs_12 = {
                'num':[0. for i in range(subtypes)],
                'length':[0. for i in range(1,subtypes)],
                'ratio':0
                }
        outputs_349 = {
                'num':[0. for i in range(subtypes)],
                'length':[0. for i in range(1,subtypes)],
                'ratio':0
                }
        outputs_50 = {
                'num':[0. for i in range(subtypes)],
                'length':[0. for i in range(1,subtypes)],
                'ratio':0,
                'kde_main':[0. for i in range(50)],
                'kde_sub':[0. for i in range(50)],
                'kde_subsub':[0. for i in range(50)]
                }
                
        outpouts_mean = {
            'year'+str(i + 1): copy.deepcopy(outputs_12) if i < 2 else copy.deepcopy(outputs_349) if i < 49 else copy.deepcopy(outputs_50)
            for i in range(N)
        }
        
        
        
        soilSpace = pb.SDF_PlantContainer(1e6, 1e6,1e6, True)  # to avoid root growing aboveground
        plant = PlantPython()

        # Open plant and root parameter from a file
        plant.readParameters( CPlantBox_dir + '/experimental/wine_fichtl/results/xmlFiles/' + genotype + "-wineV2.xml")

        plant.setGeometry(soilSpace) 


        yr_to_BEDD = 1225
        no_thin = 500
        thin_root_P = [[0.],[0.05],[1.]]
        thin_root_Page = [0., 1225 * 5, 1225 * 10]
        ps = plant.getOrganRandomParameter(1)[0]
        ps.seedPos = pb.Vector3d(0.,0.,-10.)
        for ii, pp in enumerate(plant.getOrganRandomParameter(2)):            
            if ii == 0:      
                pp.ldelay  = 0*yr_to_BEDD
                pp.ldelays = yr_to_BEDD * 50.#params['ldelays0'] #200#*5
                pp.successorP = [[[1.]]]# [[params['successorP0']]]
                pp.successorNo = [10]#[params['successorNo0']] 
                pp.successorST = [[[2]]]
                pp.successorOT = [[[2]]]
            elif ii == 1:      
                pp.ldelay  = 0*yr_to_BEDD
                pp.ldelays = yr_to_BEDD *0.5 #yr_to_BEDD * params['ldelays0'] #200#*5
                pp.successorP = [[[1]]]#[[params['successorP0']]]
                pp.successorNo = [1] #[params['successorNo0']] 
                pp.successorST = [[[2]]]
                pp.successorOT = [[[2]]]
                pp.k_survive = 1.261
                pp.lambda_survive = 9.426
                pp.rlt_winter_max = 26.8 
                pp.rlt_winter_min =  3.169
            elif ii == 2:      
                pp.ldelay_v  = [0,0]
                pp.ldelays_v = [yr_to_BEDD , yr_to_BEDD*50] 
                pp.successorNo = [1,no_thin] #[params['successorNo0']] 
                pp.successorP = [ [[0.5],[0.5],[0.1]], thin_root_P] #, [0.1, 0]]#[[params['successorP0']]]
                pp.successorST = [ [[ii + 1],[ii + 1],[ii + 1]], [[ii + 4],[ii + 4],[ii + 4]]] #[[ii + 1, ii + 4], [ii + 1, ii + 4]]
                pp.successorOT = [ [[2],[2],[2]], [[2],[2],[2]]] 
                pp.successorP_age = [[0., 1225 * 2, 1225 * 100],thin_root_Page]
                pp.k_survive = 1.261
                pp.lambda_survive = 9.426
                pp.rlt_winter_max = 26.8 
                pp.rlt_winter_min =  3.169
                pp.tropismN = 0.5  
            elif ii == 3:      
                pp.ldelay_v  = [0,0]
                pp.ldelays_v = [yr_to_BEDD , yr_to_BEDD*50]  #yr_to_BEDD * params['ldelays0'] #200#*5
                pp.successorNo = [1,no_thin] #[params['successorNo0']] 
                pp.successorP = [ [[0.5],[0.5],[0.5]],thin_root_P] # [[0.5, 0.], [0.05, 0]]#[[params['successorP0']]]
                pp.successorST = [ [[ii + 1],[ii + 1],[ii + 1]], [[ii + 4],[ii + 4],[ii + 4]]] #
                pp.successorOT = [ [[2],[2],[2]], [[2],[2],[2]]] 
                pp.successorP_age = [[0., 1225 * 2, 1225 * 10],thin_root_Page]
                # pp.r *= 2
                pp.k_survive = 1.261
                pp.lambda_survive = 9.426
                #pp.k_survive = 1.367 
                #pp.lambda_survive =  4.995 
                pp.rlt_winter_max = 13
                pp.rlt_winter_min =  0
            elif ii == 4:      
                pp.ldelay_v  = [0,0]
                pp.ldelays_v = [yr_to_BEDD , yr_to_BEDD*50] #yr_to_BEDD * params['ldelays0'] #200#*5
                pp.successorNo = [1,no_thin] #[params['successorNo0']] 
                pp.successorP = [ [[0.05]],thin_root_P] #[[0.05, 0.], [0.05, 0.]]#[[params['successorP0']]]
                pp.successorST =  [ [[ii + 1],[ii + 1]], [[ii + 4],[ii + 4],[ii + 4]]] #[[ii + 1, ii + 4], [ii + 1, ii + 4]]
                pp.successorOT = [ [[2],[2]], [[2],[2],[2]]] 
                pp.successorP_age = [[],thin_root_Page]
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
                pp.ldelays =yr_to_BEDD * 50 #yr_to_BEDD * params['ldelays0'] #200#*5
                pp.successorP = [thin_root_P]#[[params['successorP0']]]
                pp.successorNo = [no_thin] #[params['successorNo0']] 
                pp.successorST = [[[ii + 4],[ii + 4],[ii + 4]]]
                pp.successorOT = [[[2],[2],[2]]]
                pp.successorP_age = [thin_root_Page]
                #pp.k_survive = 1.261
                #pp.lambda_survive = 1.367
                #pp.rlt_winter_max = 13.9 
                #pp.rlt_winter_min =  0.
                pp.k_survive = 1.261
                pp.lambda_survive = 9.426
                pp.rlt_winter_max = 13/2/2
                pp.rlt_winter_min =  0
                
            if ii <= 5:   
                pp.a = 0.093
                pp.a_gr = 0.083/yr_to_BEDD
            else: # fine roots
                pp.r = 1. # don t know the value, i just want something high
                
                # was not used for r300
                pp.ldelay  = 0*yr_to_BEDD
                pp.ldelays =yr_to_BEDD  #yr_to_BEDD * params['ldelays0'] #200#*5
                pp.successorP = [[[1.]]]#[[params['successorP0']]]
                pp.successorNo = [1]#no_thin] #[params['successorNo0']] 
            
            
        '''
        start simulation
        '''
        # select randomly one of the rsml files of year 1
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
        #plant.betaN = 5000
          

        all_real_lengths = []
        all_real_subtypes = []
        
        all_lengths = []
        all_ages = []
        all_subtypes = []
        all_alive = []
        # all_alive_long = []
        # all_alive_short = []
        dt = yr_to_BEDD  # ~1 yr

        if doVTP:
            get3Dshape(plant,title_ = "./results/part1/vtp/"+genotype+"0", saveOnly = True)
        for i in range(N):
            print('age', i, end=", ", flush = True)
            plant.survivalTest()
            
            
            '''
            Ratio
            ''' 
            orgs_all = plant.getOrgans(2, False)
            all_real_subtypes_temp =np.array([org.param().subType for org in orgs_all  if  org.isAlive()])
            all_real_lengths_temp =np.array([org.getLength() for org in orgs_all  if  org.isAlive()])
            len_fine = 0
            is_fine_roots = np.isin(all_real_subtypes_temp,fine_root_types)
            if sum(is_fine_roots) > 0:
                len_fine = sum(all_real_lengths_temp[is_fine_roots])
            len_long = sum(all_real_lengths_temp[np.invert(is_fine_roots)]) # ignore the length of 1 as seem different between xml and xlsx
            value = len_fine/(len_long + len_fine)
            print('percent_after_winter', np.round(value*100), end=", ", flush = True)   
            
            plant.simulate(dt, False) # i * 
            
            orgs_all = plant.getOrgans(2, False)
            orgs_ = []
            for oo in orgs_all:
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
            
            
            all_real_subtypes.append(np.array([org.param().subType for org in orgs_all  if  org.isAlive()]))
            all_real_lengths.append(np.array([org.getLength() for org in orgs_all  if  org.isAlive()]))
            
            '''
            Ratio
            ''' 
            len_fine = 0
            is_fine_roots = np.isin(all_real_subtypes[i],fine_root_types)
            if sum(is_fine_roots) > 0:
                len_fine = sum(all_real_lengths[i][is_fine_roots])
            len_long = sum(all_real_lengths[i][np.invert(is_fine_roots)]) # ignore the length of 1 as seem different between xml and xlsx
            value = len_fine/(len_long + len_fine)
            print('percent', np.round(value*100), end=", ", flush = True)    
            '''
            Num
            '''
            #nums = []
            #for st in range(subtypes): # from 2 to 5            
            #    nums.append(sum(all_alive[i][all_subtypes[i] == (st + 1)]))
            #nums.append(sum(all_alive[i][all_subtypes[i] > (st + 1)]))
            # print('num', nums, end=", ", flush = True)
            # all_alive_long.append(np.array([org.isAlive() for org in orgs_ if org.param().subType in long_root_types]))
            # all_alive_short.append(np.array([org.isAlive() for org in orgs_ if org.param().subType in fine_root_types]))

            # # get3Dshape(plant,title_ = 'wine'+str(i+1), saveOnly = True) 
            if doVTP:
                get3Dshape(plant,title_ = "./results/part1/vtp/"+genotype+str(i+1), saveOnly = True)
            
        # print('postprocessing')
        #print("--- %s seconds for plant development---" % (time.time() - start_time),rep)
        
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
            #outpouts_mean['year'+str(year+1)]['num'][subtypes] = sum(all_alive[year][all_subtypes[year] > subtypes])
            '''
            Ratio
            ''' 
            len_fine = 0
            is_fine_roots = np.isin(all_real_subtypes[year],fine_root_types)
            if sum(is_fine_roots) > 0:
                len_fine = sum(all_real_lengths[year][is_fine_roots])
            len_long = sum(all_real_lengths[year][np.invert(is_fine_roots)]) # ignore the length of 1 as seem different between xml and xlsx
            value = len_fine/(len_long + len_fine)
            outpouts_mean['year'+str(year+1)]['ratio'] = value     
            

        '''
        Age distribution at last time step
        '''
        idy = N - 1
        
             
        for key, value in orders.items():
            # List to store KDE y-values for each FileName
            
            x = xx[key]
            
            root_age = all_ages[np.isin(all_subtypes[idy],value)]
            try:
                if len(root_age) == 0:
                    y = np.zeros(len(x))
                elif len(root_age) == 1:
                    y = np.zeros(len(x))
                    indx_y = np.where(abs(x - root_age[0]) == min(abs(x - root_age[0])))[0][0]
                    y[indx_y] = 1.
                else:
                    kde = gaussian_kde(root_age)
                    y = kde(x)
            except:
                print('issue root_age',root_age)
                raise Exception
            
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
    
        print('output',output)
        print("--- %s seconds for plant development---" % ((time.time() - start_time)/(rep+1)))
        # print("--- %s seconds for total---" % (time.time() - start_time),'at',(rep+1),'/',(reps))
        
    #if getdict:
    #    return outpouts_mean
        
    # print('outpouts_mean',outpouts_mean)
    #f = [i for s in output for i in s]
    return output #np.array(f)
        


    
if __name__ == '__main__': 
    
    genotype = sys.argv[1]
    rep = int(sys.argv[2])
    
    extraName = "default"
    if len(sys.argv) > 3:
        extraName = sys.argv[3]
    print('sys.argv',sys.argv)
    if rep == 0:
        directory ='./results/outputSim/'+extraName
        os.makedirs(directory, exist_ok=True)
        
    with open(CPlantBox_dir + '/experimental/wine_fichtl/results/objectiveData/measurements'+ genotype +'InitXX.pkl','rb') as f:
        xx = pickle.load(f)
        
    output = run_benchmark(xx, genotype, rep)
    
    with open(CPlantBox_dir + '/experimental/wine_fichtl/results/outputSim/'+extraName+'/'+ genotype + str(rep) +'.pkl','wb') as f:
         pickle.dump(output,f, protocol=pickle.HIGHEST_PROTOCOL)
            
    #comm.Barrier()
    #if rank == 0:
    #    import gatherAndPlotOutputs
    #    gatherAndPlotOutputs.plot_all()
    
