'''
Create the objective file for each genotype.
[mean num subtype 1 year 1]
..
[mean num subtype 5 year 1]
..
[mean num subtype 1 year 10]
..
[mean num subtype 5 year 10]
[mean sum length subtype 1 year 2] !! not length of num 1
..
[mean sum length 5 year 1]
..
[mean sum length 1 year 10]
..
[mean sum length 5 year 10]
[mean ratio length fine roots year 1] (ratio => except from length of 1)
..
[mean ratio length fine roots year 10]
[sd num subtype 1 year 1]
..
[sd num subtype 5 year 1]
..
[sd num subtype 1 year 10]
..
[sd num subtype 5 year 10]
[sd sum length subtype 1 year 1]
..
[sd sum length 5 year 1]
..
[sd sum length 1 year 10]
..
[sd sum length 5 year 10]
[sd ratio length fine roots year 1]
..
[sd ratio length fine roots year 10]
'''
import sys; sys.path.append("../.."); sys.path.append("../../src/")
import plantbox as pb
import pickle

import matplotlib.pyplot as plt
import numpy as np
import os
from structural.Plant import PlantPython
import pandas as pd
import copy


def getMeasData(genotype,dt_types, years_ = 10):#, doGraphs = False):
    #genotype = "B" # "D" "E"
    subtypes = 5


    # Define column names
    columns = ['Year', 'Data', 'ST']

    # Create an empty DataFrame with these columns
    dfmean = pd.DataFrame(columns=columns)
    dfsd = pd.DataFrame(columns=columns)

    # outpouts = [
        # [
            # [[0 for ii in range(subtypes)] for yy in range(years_)],  # first element
            # [[0 for ii in range(subtypes)] for yy in range(years_)],  # second element
            # [0 for yy in range(years_)]                               # third element
        # ]
        # for sdmn in range(2)
    # ]
    # # [[5 [subtype] * 10 [yrs] * 2 [data type],[10 [yrs] * 1 [data type]]] * 2 [mean + sd]]

    outputs_12 = {
            'num':[0.,0.,0.,0.,0.],
            'length':[0.,0.,0.,0.],
            'ratio':0
            }
    outputs_3to10 = {
            'num':[0.,0.,0.,0.,0.],
            'ratio':0
            }
            
    outpouts_mean = {
        'year'+str(i+1): copy.deepcopy(outputs_12) if i < 2 else copy.deepcopy(outputs_3to10)
        for i in range(10)
    }
    outpouts_sd = {
        'year'+str(i+1): copy.deepcopy(outputs_12) if i < 2 else copy.deepcopy(outputs_3to10)
        for i in range(10)
    }

    for year in [1,2]:

        file_path = "./rsml/RSML_year"+str(year)+"/"
        file_names = [name for name in os.listdir(file_path) if name[0] == genotype ] # repetitions of same genotype
        
        """
        RSMLs:
        Type 1: the initially planted roots
        Type 2: first order roots (begin to grow at once)
        Type 3: second order roots (emergence times are estimated based on parametrisation of first order roots) 
        """
        #outpouts_ = np.zeros((len(file_names), 3, subtypes)) #np.array([ [ [0 for ii in range(subtypes)]   ,[0 for ii in range(subtypes)]  , [0] ]  for i in range(len(file_names))])
        outpouts__ = {
            'file'+str(i+1): copy.deepcopy(outputs_12) for i in range(len(file_names))
            }
        for nm_id, name in enumerate(file_names):
            print(name)
            name = name[:-5]
            plant = PlantPython()
            path_ = "../../modelparameter/structural/rootsystem/"
            plant.readParameters(path_ +  "wine_Fichtl.xml")
            radii, types, lengths = plant.get_rsml_data(file_path + name + ".rsml") 
            #print('lengths', lengths[types == 1], sum(lengths[types == 1]))
            
            # is_fine_roots = np.zeros(types.shape)
            
            
            # and also do not have long lived roots attached to them?
            if year == 1:
                is_fine_roots = (lengths < 5) & (radii < 0.093/2.)# & (types > 1) # (lengths < 6) & because we have extra restrictions, more roots go into that category. better set an 'or' maybe
            else:
                is_fine_roots = (lengths < 5)# & (types > 1)
            
            types = np.array([st if is_fine_roots[i] == 0 else -1 for i, st in enumerate(types)])   


            #print('sum fine roots',sum(is_fine_roots), lengths, radii)
            for dt_id, datatype in enumerate(['num', 'length', 'ratio']):
                #print('dt_id', nm_id, dt_id, len(outpouts__), len(outpouts__['file'+str(nm_id+1)]), outpouts__['file'+str(nm_id+1)][datatype])
                if datatype == 'num':
                    for st in range(subtypes):
                        value = sum(types == (st + 1))
                        #outpouts[0][year-1][st] += num_st/len(file_names)
                        #print('outpouts_value', st,outpouts_[nm_id][dt_id][st], value, 'get',outpouts_)
                        # outpouts_[nm_id][dt_id][st] = value
                        outpouts__['file'+str(nm_id+1)][datatype][st] = value
                        #print('outpouts_value', st,outpouts_[nm_id][dt_id][st], value, 'get',outpouts_)
                        
                elif datatype == 'length':
                    for st in range(1,subtypes):
                        value = sum(lengths[types == (st + 1)])
                        #outpouts[0][year-1][st] += sum_len/len(file_names)
                        #outpouts_[nm_id][dt_id][st] = value
                        outpouts__['file'+str(nm_id+1)][datatype][st-1] = value
                else:
                    if sum(is_fine_roots) > 0:
                        len_fine = sum(lengths[is_fine_roots])
                    else:
                        len_fine = 0
                    len_long = sum(lengths[types > 1]) # ignore the length of 1 as seem different between xml and xlsx
                    value = len_fine/(len_long + len_fine)
                    #outpouts_[nm_id][dt_id][0] = value
                    outpouts__['file'+str(nm_id+1)][datatype] = value
                
                    
        #print('outpouts_',outpouts_)        
        for dt_id, datatype in enumerate(['num', 'length', 'ratio']):
            if datatype == 'ratio':
                values = [outpouts__['file'+str(nm_id+1)][datatype] for nm_id, name in enumerate(file_names)]
                outpouts_mean['year'+str(year)][datatype] = sum(values)/len(file_names)
                # outpouts[0][dt_id][year-1] = sum(outpouts_[:,dt_id,0])/len(file_names)
                outpouts_sd['year'+str(year)][datatype] = np.sqrt(sum((values - outpouts_mean['year'+str(year)][datatype]) ** 2)/len(file_names))
            elif datatype == 'length':
                for st in range(1,subtypes):
                    values = np.array([outpouts__['file'+str(nm_id+1)][datatype][st-1] for nm_id, name in enumerate(file_names)])
                    
                    outpouts_mean['year'+str(year)][datatype][st-1] = sum(values)/len(file_names)
                    outpouts_sd['year'+str(year)][datatype][st-1] = np.sqrt(sum((values - outpouts_mean['year'+str(year)][datatype][st-1]) ** 2)/len(file_names))
            else:
                for st in range(subtypes):
                    values = np.array([outpouts__['file'+str(nm_id+1)][datatype][st] for nm_id, name in enumerate(file_names)])
                    outpouts_mean['year'+str(year)][datatype][st] = sum(values)/len(file_names)
                    outpouts_sd['year'+str(year)][datatype][st] = np.sqrt(sum((values - outpouts_mean['year'+str(year)][datatype][st]) ** 2)/len(file_names))
    
    
    #print('means',outpouts[0])
    #print('sd',outpouts[1])

    type1_yr1 = outpouts_mean['year1']['num'][0] # outpouts[0][0][0][0]
    type1_yr2 = outpouts_mean['year2']['num'][0] #outpouts[0][0][1][0]
    dt_type1 = type1_yr2 - type1_yr1
    if dt_type1 >= 0:
        dt_type1 = np.mean(dt_types)

    # # length of type 1 will not be used
    # outpouts[0][0][0][0] = 0.
    # outpouts[1][0][0][0] = 0.
    # outpouts[0][1][0][0] = 0.
    # outpouts[1][1][0][0] = 0.


    ratio_init = outpouts_mean['year1']['ratio'] #outpouts[0][2][0]
    sd_init = outpouts_sd['year1']['ratio']
    # CV_all = sd_init/ratio_init # makes it too big by the end

    #print([outpouts_mean['year'+str(yr)]['ratio'] for yr in [1,2]])
    #raise Exception
    mean_ratio_init = 0.08
    ratio_10 = 0.8 - mean_ratio_init + ratio_init

    year = 2
    dt_id = 2
    # outpouts[0][dt_id][year-1] = ratio_init + (ratio_10 - ratio_init ) * (year/10)
    # outpouts[1][dt_id][year-1] = CV_all * outpouts[0][dt_id][year-1] 
    outpouts_mean['year'+str(year)]['ratio'] = ratio_init + (ratio_10 - ratio_init ) * ((year-1)/9)
    outpouts_sd['year'+str(year)]['ratio'] = sd_init #CV_all * outpouts_mean['year'+str(year)]['ratio']

    for year in range(3,11):
        for dt_id, datatype in enumerate(['num', 'ratio']): # 'length',
            if datatype == 'ratio':#'num':
                # outpouts[0][dt_id][year-1] = ratio_init + (ratio_10 - ratio_init ) * (year/10)
                # outpouts[1][dt_id][year-1] = CV_all * outpouts[0][dt_id][year-1] 
                outpouts_mean['year'+str(year)][datatype] = ratio_init + (ratio_10 - ratio_init ) * ((year-1)/9)
                outpouts_sd['year'+str(year)][datatype] = sd_init# CV_all * outpouts_mean['year'+str(year)][datatype]
            elif datatype == 'length':
                raise Exception
                for st in range(1,subtypes):
                    outpouts_mean['year'+str(year)][datatype][st-1] = outpouts_mean['year2'][datatype][st-1]
                    outpouts_sd['year'+str(year)][datatype][st-1] = outpouts_sd['year2'][datatype][st-1]
            else:
                for st in range(subtypes):
                    #outpouts[0][dt_id][year-1][st] = outpouts[0][dt_id][1][st]
                    #outpouts[1][dt_id][year-1][st] = outpouts[1][dt_id][1][st]
                    outpouts_mean['year'+str(year)][datatype][st] = outpouts_mean['year2'][datatype][st]
                    outpouts_sd['year'+str(year)][datatype][st] = outpouts_sd['year2'][datatype][st]
                            
                # outpouts[0][0][year-1][0] = max(0.,outpouts[0][0][year-1][0] + dt_type1 * year)
                outpouts_mean['year'+str(year)][datatype][0] = max(0.,outpouts_mean['year'+str(year)][datatype][0] + dt_type1 * float(year - 2.))
     
    # for dt_id, datatype in enumerate(['num', 'length', 'ratio']):
        # print(datatype)
        # print([outpouts_mean['year'+str(year)][datatype] for year in range(1,11)])
    outpouts = {'mean':outpouts_mean, 'sd':outpouts_sd}
    
    with open('./measurements'+ genotype +'Init.pkl','wb') as f:
         pickle.dump(outpouts,f, protocol=pickle.HIGHEST_PROTOCOL)
         
    def flatten_values(obj):
        if isinstance(obj, dict):
            for v in obj.values():
                yield from flatten_values(v)
        elif isinstance(obj, (list, tuple, np.ndarray)):
            for v in obj:
                yield from flatten_values(v)
        else:
            yield obj
    mu_tmp = np.array(list(flatten_values(outpouts_mean)))        
    #mu_tmp = np.concatenate([np.ravel(list(oo.values())) for oo in outpouts_mean.values()])
    sd_tmp = np.array(list(flatten_values(outpouts_sd)))   #np.ravel([list(oo.values()) for oo in outpouts_sd.values()])
    print('mu_tmp',mu_tmp, len(mu_tmp))
    # first random guess
    sd_mu = mu_tmp * 0.1  
    sd_sd = sd_tmp * 0.1    
    #sd >= 0.01 to avoid division by 0 when computing the log.
    #sd_mu[sd_mu <0.01] = 0.01
    # sd_sd[sd_sd < 0.01] = 0.01

    with open('./measurements'+ genotype +'.pkl','wb') as f:
         pickle.dump({'sd_sd':sd_sd,
         'sd_data':sd_tmp,
         'mu_data' :mu_tmp,
         'sd_mu' : sd_mu},f, protocol=pickle.HIGHEST_PROTOCOL)

    
    return dt_type1 #outpouts_mean, outpouts_sd

dt_types = []
for genotype in ["B", "D", "E"]:
    dt_types.append(getMeasData(genotype, dt_types))
print(dt_types)

