
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
    for scenario in [0,1]:        

        with open('./results/outputSim/'+extra_name+'/'+ genotype    + '_soil_' + str(scenario) +'.pkl','rb') as f:
            output_all = pickle.load(f)

        # scaling = np.unique(np.logspace(np.log10(1), np.log10(50), num=10).astype(int)) 
        scaling = np.array([ 1,  2,  3,  5,  8, 13, 20, 32, 50])
        output = [out for year, out in enumerate(output_all)]# if year in scaling]
        times = [out['time'] for out in output]
        transpirations = [out['transpiration'] for out in output]
        h_bs = [out['h_bs'] for out in output]
        qs = [out['q'] for out in output]
        depths = [out['depth_array'] for out in output]

        print(len(output), len(scaling))
        if True:
            fig, ax1 = plt.subplots()
            dt_soil = 3600. / (24 * 3600)
            for id_, y_ in enumerate(transpirations): 
                ax1.plot(times[id_],-np.array(y_),  label = 'year ' + str(scaling[id_]))  # cumulative transpiratio
            ax1.set_xlabel("Time [d]")
            ax1.legend()
            fig.suptitle("Transpiration $[mL/day]$ per plant", fontsize=16)
            ax1.grid() 
            plt.tight_layout()
            plt.savefig("./results/part1/"+genotype+"/transpiration"+ extra_name + str(scenario)+".jpg")
            plt.close()
            print("did Transpiration")

            fig, ax1 = plt.subplots()
            dt_soil = 3600. / (24 * 3600)
            for id_, y_ in enumerate(transpirations): 
                ax1.plot(times[id_], np.cumsum(-np.array(y_) * dt_soil),  label = 'year ' + str(scaling[id_]))  # cumulative transpiratio
            ax1.set_xlabel("Time [d]")
            ax1.legend()
            fig.suptitle("Cumulative Transpiration $[mL]$ per plant", fontsize=16)
            ax1.grid() 
            plt.tight_layout()
            plt.savefig("./results/part1/"+genotype+"/Cumultranspiration"+ extra_name + str(scenario)+".jpg")
            plt.close()
            print("did cumulative Transpiration")

        fig, ax = plt.subplots(1,1, figsize=(5,10))
        for id_, h_bs_ in enumerate(h_bs): 
            ax.plot(h_bs_, depths[id_],  label = 'year ' + str(scaling[id_]))
        ax.legend()
        ax.grid() 
        fig.suptitle('Soil water potential [cm]', fontsize=16)
        plt.tight_layout()
        plt.savefig("./results/part1/"+genotype+"/h_bs"+ extra_name + str(scenario)+".jpg")
        plt.close()
        print("did h_bs fraction")
        
if __name__ == "__main__":

    extra_name = "defaultsoil"
    if len(sys.argv) > 1:
        extra_name = sys.argv[1]
    print("using folder", extra_name)
    scaling = np.array([ 1,  2,  3,  5,  8, 13, 20, 32, 50])
    if True:
        for genotype in ['B', 'D', 'E']:
            #directory = "./results/part1/"+genotype
            #os.makedirs(directory, exist_ok=True)

            for scenario in [0,1]:
                output = []
                for i in scaling:#range(50):
                    try:
                        print('year', i)
                        
                        with open( './results/outputSim/'+extra_name+'/'+ genotype + '0_soil_' + str(scenario) + "_"  + str(i)+'.pkl','rb') as f:
                            temp = pickle.load(f)
                            output.append(temp)
                    except:
                        print( './results/outputSim/'+extra_name+'/'+ genotype + '0_soil_' + str(scenario) + "_"  + str(i))
                        print("skip",genotype, scenario, i)
                        raise Exception

                with open(  './results/outputSim/'+extra_name+'/'+ genotype  + '_soil_' + str(scenario) + '.pkl','wb') as f:
                    #output = [item for sublist in output for item in sublist]
                    pickle.dump(output,f, protocol=pickle.HIGHEST_PROTOCOL)
    if True:
        print('extra_name is',extra_name)
        print("\tdo B")
        compaireOutPuts('B', extra_name)
        print("\n\n\tdo D")
        compaireOutPuts('D', extra_name)
        print("\n\n\tdo E")
        compaireOutPuts('E', extra_name)

    