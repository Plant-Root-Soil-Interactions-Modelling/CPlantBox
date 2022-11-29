""" water movement within the root (static soil) """


import sys; 
import math
import os
import numpy as np
import pandas as pd

def toTry():
    
    dfP2 = pd.read_csv('results/forP2.csv')


    # dictPara= {'Qsv' : dfP2['Qsv'],
    #                 'nodeDv':dfP2['nodeDv'], 
    #            'MulimSucv':dfP2['MulimSucv'],
    #                 'kaav':dfP2['kaav'], 
    #            'kssv':dfP2['kssv'],
    #           'GrRatiov':dfP2['GrRatiov'],
    #           'CarbonCostv':dfP2['CarbonCostv']}
    # df = pd.DataFrame(data=dictPara)
    if sum(dfP2.duplicated()) > 0 :
        print("CalibP2Database::toTry(): sum(df.duplicated()) > 0")
        raise Exception
    return dfP2 #dictPara

