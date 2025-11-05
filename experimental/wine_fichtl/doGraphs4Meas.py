import sys; sys.path.append("../.."); sys.path.append("../../src/")
import plantbox as pb
import pickle

import matplotlib.pyplot as plt
import numpy as np
import os
from structural.Plant import PlantPython
import pandas as pd
import copy

###
# Create the graphs
# to do: test the survival graphs + update code. as long as attached to existng long lived is ok
# as long as a living root is attched, will keep on living.
##
outpout_mean = {}
outpout_sd = {}
for genotype in ["B", "D", "E"]:

    with open('./measurements'+ genotype +'Init.pkl','rb') as f:
        temp = pickle.load(f)
        outpout_mean[genotype] = temp['mean']
        outpout_sd[genotype] = temp['sd']
        

#####
#   num
#####
# create dict for each genotype, 2D arrays: [[length for each year] each subtype]
outpout_num_mean = {}
outpout_num_sd = {}

for genotype in ["B", "D", "E"]:
    outpout_num_mean[genotype] = np.array([
            np.array([outpout_mean[genotype]['year'+str(year + 1)]['num'][st] for year in range(10)])
            for st in range(subtypes)])
            
        