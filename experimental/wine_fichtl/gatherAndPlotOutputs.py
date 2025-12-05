
import sys; 
from pathlib import Path

#edit CPlantBox dir
CPlantBox_dir =  "../.." # "/data2model_0807/CPlantBox"

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


subtypes_long = 5
subtypes = 9
years = 50
ages = [i+1 for i in range(50)]

    
def compaireOutPuts(genotype, extra_name):
            
    with open('./results/objectiveData/measurements'+ genotype +'InitXX.pkl','rb') as f:
        xx = pickle.load(f)
            
    with open('./results/objectiveData/measurements'+ genotype +'Init.pkl','rb') as f:
        temp = pickle.load(f)
        obs_mean = temp['mean']
        obs_sd = temp['sd']
        
    with open('./results/outputSim/'+extra_name+'/'+ genotype   +'.pkl','rb') as f:
        output = pickle.load(f)
        nums = [[[out['year'+str(year+1)]['num'][st] for year in range(50)] for st in range(subtypes)] for out in output]
    
    fig, axs = plt.subplots(2, 2, figsize=(12, 8))
    for st_, ax in enumerate(axs.flat):
        st = st_ + 5
        for _, num in enumerate(nums):
            if _ == 0:
                ax.plot(ages, num[st], color='black', label = 'sim')
            else:
                ax.plot(ages, num[st], color='black')
        ax.legend()
        ax.grid()
        ax.set_title(st+1)
    fig.suptitle('Number of THIN roots per order--50yrs', fontsize=16)
    plt.tight_layout()
    plt.savefig("./results/part1/"+genotype+"/THINNum50"+ extra_name +".jpg")
    plt.close()
    print("did Number of THIN roots per order--50yrs")
    
    ratios = [[out['year'+str(year+1)]['ratio'] for year in range(50)] for out in output]
    fig, ax = plt.subplots(1, 1, figsize=(10, 4))
    for _, ratio in enumerate(ratios):
        if _ == 0:
            ax.plot(ages, ratio,  color='black', label = 'sim')
        else:
            ax.plot(ages, ratio,  color='black')
    ax.plot([1,10,50], [0.05,0.8,0.8],  color='blue', label = 'objective')
    ax.legend()
    ax.grid()
    fig.suptitle('fine:total root length ratio', fontsize=16)
    plt.savefig("./results/part1/"+genotype+"/Ratio"+ extra_name +".jpg")
    plt.close()
    print("did Ratio")
    
    fig, axs = plt.subplots(2, 3, figsize=(12, 8))
    for st, ax in enumerate(axs.flat):
        if st < 5:
            obs_age = [1, 2, 50]
            y = np.array([obs_mean['year' + str(y)]['num'][st] for y in obs_age])
            error = np.array([obs_sd['year' + str(y)]['num'][st] for y in obs_age]) #obs_sd['num'][st]
            #axs[gid].plot(x, y, label = st + 1)
            ax.fill_between(obs_age, y - error, y + error, color='blue', alpha=0.1, label = 'obs')
            for _, num in enumerate(nums):
                if _ == 0:
                    ax.plot(ages, num[st], color='black', label = 'sim')
                else:
                    ax.plot(ages, num[st], color='black')
            ax.legend()
            ax.grid()
            ax.set_xlim(right=2.5)
            ax.set_title(st+1)
        else:
            ax.axis('off')
    fig.suptitle('Number of roots per order--2yrs', fontsize=16)
    plt.tight_layout()
    plt.savefig("./results/part1/"+genotype+"/Num2"+ extra_name +".jpg")
    plt.close()
    print("did Number of roots per order--2yrs")
    
    fig, axs = plt.subplots(2, 3, figsize=(12, 8))
    for st, ax in enumerate(axs.flat):
        if st < 5:
            obs_age = [1, 2, 50]
            y = np.array([obs_mean['year' + str(y)]['num'][st] for y in obs_age])
            error = np.array([obs_sd['year' + str(y)]['num'][st] for y in obs_age]) #obs_sd['num'][st]
            #axs[gid].plot(x, y, label = st + 1)
            ax.fill_between(obs_age, y - error, y + error, color='blue', alpha=0.1, label = 'obs')
            for _, num in enumerate(nums):
                if _ == 0:
                    ax.plot(ages, num[st], color='black', label = 'sim')
                else:
                    ax.plot(ages, num[st], color='black')
            ax.legend()
            ax.grid()
            ax.set_title(st+1)
        else:
            ax.axis('off')
    fig.suptitle('Number of roots per order--50yrs', fontsize=16)
    plt.tight_layout()
    plt.savefig("./results/part1/"+genotype+"/Num50"+ extra_name +".jpg")
    plt.close()
    print("did Number of roots per order--50yrs")
    
    
    
    lengths = [[[out['year'+str(year+1)]['length'][st-1] for year in range(50)] for st in range(1,subtypes)] for out in output]
    
    fig, axs = plt.subplots(2, 2, figsize=(8, 8))
    for st, ax in enumerate(axs.flat):
        obs_age = [1, 2]
        y = np.array([obs_mean['year' + str(y)]['length'][st] for y in obs_age])
        error =  np.array([obs_sd['year' + str(y)]['length'][st] for y in obs_age])
        #axs[gid].plot(x, y, label = st + 1)
        ax.fill_between([1, 2], y - error, y + error, color='blue', alpha=0.1, label = 'obs')
        for _, length in enumerate(lengths):
            if _ == 0:
                ax.plot(ages, length[st], color='black', label = 'sim')
            else:
                ax.plot(ages, length[st], color='black')
        ax.legend()
        ax.set_xlim(right=2.5)
        ax.grid()
        ax.set_title(st+2)
    fig.suptitle('Total root length per order--2yrs', fontsize=16)
    plt.tight_layout()
    plt.savefig("./results/part1/"+genotype+"/Len2"+ extra_name +".jpg")
    plt.close()
    print("did Total root length per order--2yrs")
    
    fig, axs = plt.subplots(2, 2, figsize=(8, 8))
    for st, ax in enumerate(axs.flat):
        obs_age = [1, 2]
        y = np.array([obs_mean['year' + str(y)]['length'][st] for y in obs_age])
        error =  np.array([obs_sd['year' + str(y)]['length'][st] for y in obs_age])
        #axs[gid].plot(x, y, label = st + 1)
        ax.fill_between([1, 2], y - error, y + error, color='blue', alpha=0.1, label = 'obs')
        for _, length in enumerate(lengths):
            if _ == 0:
                ax.plot(ages, length[st], color='black', label = 'sim')
            else:
                ax.plot(ages, length[st], color='black')
        ax.legend()
        ax.grid()
        ax.set_title(st+2)
    fig.suptitle('Total root length per order--50yrs', fontsize=16)
    plt.tight_layout()
    plt.savefig("./results/part1/"+genotype+"/Len50"+ extra_name +".jpg")
    plt.close()
    print("did Total root length per order--50yrs")
    
    
    fig, axs = plt.subplots(1, 3, figsize=(10, 4))
    for st, order in enumerate(['main', 'sub', 'subsub']):        
        mean_kde = obs_mean['year50']['kde_' + order]
        std_kde = obs_sd['year50']['kde_' +order] 
        # axs[st].plot(xx[order],mean_kde, label='obs', linewidth=2)
        axs[st].fill_between(xx[order], mean_kde - std_kde, mean_kde + std_kde, color='gray', alpha=0.3)
        for _, out in enumerate(output):
            if _ == 0:
                axs[st].plot(xx[order], out['year50']['kde_' + order],  color='black', label = 'sim')
            else:
                axs[st].plot(xx[order], out['year50']['kde_' + order],  color='black')
        axs[st].grid()
        axs[st].set_title(order)
    fig.suptitle('Survival probability', fontsize=16)
    plt.tight_layout()
    plt.savefig("./results/part1/"+genotype+"/Density"+ extra_name +".jpg")
    plt.close()
    print("did Survival probability")
  

if __name__ == "__main__":
    extra_name = sys.argv[1]

    for genotype in ['B', 'D', 'E']:
        directory = "./results/part1/"+genotype
        os.makedirs(directory, exist_ok=True)
        output = []
        for rep in range(8):

            with open( './results/outputSim/'+extra_name+'/'+ genotype + str(rep) +'.pkl','rb') as f:
                temp = pickle.load(f)
                output.append(temp)
            print(genotype, rep, len(output))

        with open(  './results/outputSim/'+extra_name+'/'+ genotype + '.pkl','wb') as f:
            output = [item for sublist in output for item in sublist]
            pickle.dump(output,f, protocol=pickle.HIGHEST_PROTOCOL)

    print('extra_name is',extra_name)
    print("\tdo B")
    compaireOutPuts('B', extra_name)
    print("\n\n\tdo D")
    compaireOutPuts('D', extra_name)
    print("\n\n\tdo E")
    compaireOutPuts('E', extra_name)

    