import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import pandas as pd
import pandas as pd
import numpy as np

font = {'size'   : 14}
matplotlib.rc('font', **font)


df = pd.read_csv("../../data_magda/diam_length_exudates.csv")
print(df.columns.values.tolist()) 
DAS = df["DAS"].loc[:].values

soils = ['L', 'S']
gt = ['WT', 'rth3']

#colors
a = [0, 0, 1]
b = [106/255, 166/255, 1]
c = [1, 0, 0]
d = [1, 128/255, 0]
cols = [a,b,c,d]
label_root = ["L_WT", "L_rth3", "S_WT", "S_rth3"]

label_exud = [['a', 'a','a','a'],[ 'a','ab', 'a','ab'],[ 'a','b','a','b'],[ 'a','b','a','b']]


fig, ax = plt.subplots(3, figsize = (18, 10))

for i in range(0,len(soils)):
    for j in range(0, len(gt)):

        y = df[(df['Substrate']== soils[i]) & (df['Genotype']==gt[j])]

        ax[0].plot(y['DAS'], y['Phenols'], color = cols[i+2*j], label = soils[i]+'_'+gt[j])
        ax[0].fill_between(y['DAS'], (y['Phenols']-y['Phenols std']/2), (y['Phenols']+y['Phenols std']/2), color = cols[i+2*j], alpha = 0.2)
        ax[0].set_ylabel('Phenol exudation rate \n (mmol / plant / h)')

        ax[1].plot(y['DAS'], y['Sugars'], color = cols[i+2*j], label = soils[i]+'_'+gt[j])
        ax[1].fill_between(y['DAS'], y['Sugars']-y['Sugars std']/2, y['Sugars']+y['Sugars std']/2, color = cols[i+2*j], alpha = 0.2)
        ax[1].set_ylabel('Sugar exudation rate \n (mmol / plant / h)')

        ax[2].plot(y['DAS'], y['AminoAcids'], color = cols[i+2*j], label = soils[i]+'_'+gt[j])
        ax[2].fill_between(y['DAS'], y['AminoAcids']-y['AminoAcids std']/2, y['AminoAcids']+y['AminoAcids std']/2, color = cols[i+2*j], alpha = 0.2)
        ax[2].set_ylabel('Amino acid exudation rate\n (mmol / plant / h)')

        # zip joins x and y coordinates in pairs
        #label = label_exud[i+2*j]
        #kk = 0
        #for x,y in zip(y['DAS'], y['exud']):
        #    ax[2].annotate(label[kk],(x,y),textcoords="offset points",xytext=(0,10),ha='center')
        #    kk = kk+1
        

ax[0].legend()
for i in range(0,3):
    ax[i].set_xlabel('Time (d)')

plt.tight_layout()
plt.show()



