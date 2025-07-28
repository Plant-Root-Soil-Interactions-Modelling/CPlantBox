import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import visualisation.vtk_plot as vp
import matplotlib.pyplot as plt
import numpy as np

path = "../../modelparameter/structural/rootsystem/"
name = "wheat"

months=8 #|\label{l2_3:timebegin}|
times=np.linspace(0,30*months,months+1) #|\label{l2_3:timeend}|

# 72 cm*45 cm size plot
M = 16  # number of plants in rows #|\label{l2_3:plantsetbegin}|
N = 7 # number of rows
distp = 3  # distance between the root systems along row[cm]
distr =12  # distance between the rows[cm]
interrow=(M-1)*distp # intra-row spacing
row=(N-1)*distr # row spacing #|\label{l2_3:plantsetend}|

r, depth, layers = 4.2/2, 160., 32 # Soil core analysis #|\label{l2_3:soilcorebegin}|
layerVolume = depth / layers * r * r * np.pi
z_ = np.linspace(0, -depth, layers)  # slices of soil core

soilcolumn = pb.SDF_PlantContainer(r, r, depth, False)#square = false

soilcor_x = interrow/4
soilcor_y = distr
x_= [11.25, 22.5, 33.75, 11.25, 22.5, 33.75, 11.25, 22.5, 33.75, 11.25, 22.5, 33.75, 11.25, 22.5, 33.75]                    
y_= [66.0, 66.0, 66.0, 54.0, 54.0, 54.0, 42.0, 42.0, 42.0, 30.0, 30.0, 30.0, 18.0, 18.0, 18.0] 
soilcolumns_ =[pb.Vector3d(x_ij, y_ij,0) for x_ij, y_ij in zip(x_, y_)]
soilcolumns = [pb.SDF_RotateTranslate(soilcolumn, vi) for vi in soilcolumns_] #|\label{l2_3:soilcoreend}|

soilSpace = pb.SDF_PlantContainer(500,500, 500, True)

allRS = []
for i in range(0, M): #|\label{l2_3:simulationbegin}|
    for j in range(0, N):
        plant = pb.Plant()
        plant.readParameters(path + name + ".xml", fromFile = True, verbose = False )
        seed = plant.getOrganRandomParameter(pb.seed)[0]
        seed.seedPos = pb.Vector3d(distp * i, distr * j, -3.)  # cm
        plant.setGeometry(soilSpace)
        plant.initialize(False)
        plant.simulate(30*months, False)
        allRS.append(plant) 
        if i+ j == 0:
            allAna = pb.SegmentAnalyser(plant)
        else:
            allAna.addSegments(plant) #|\label{l2_3:simulationend}|
rld=np.zeros([len(soilcolumns)*len(times[1:]), layers])

for k in range(len(soilcolumns)): #|\label{l2_3:soilcolselectbegin}|
    ana = pb.SegmentAnalyser(allAna)  
    ana.crop(soilcolumns[k]) #select soil column
    for j in range(len(times[1:])):
        ana.filter("creationTime", 0, np.flip(np.asarray(times))[j]) 
        ana.pack()
        distrib = ana.distribution("length", 0., -depth, layers, True)
        rld[(len(times[1:])-1-j)*len(soilcolumns)+k]= np.array(distrib)/layerVolume #|\label{l2_3:soilcolselectend}|
    
pt_idx=[(0,0),(0,1),(0,2), #|\label{l2_3:plotbegin}|
(1,0),(1,1),(1,2),
(2,0),(2,1),(2,2),
(3,0),(3,1),(3,2),
(4,0),(4,1),(4,2),]
legend_lst = [str(int(t_i)) for t_i in times[1:]]    
fig, axes = plt.subplots(nrows = 5, ncols = int(len(soilcolumns)/5), sharex=True, sharey=True,figsize = (8,16))

for k in range(len(soilcolumns)):
    axes[pt_idx[k]].set_title('Soil core'+' ' +str(k+1))
    for j in range(len(times[1:])):
        axes[pt_idx[k]].plot(np.array(rld[len(soilcolumns)*j+k]), z_)
        axes[pt_idx[k]].set_xlim(0,5)
        
plt.setp(axes[-1, :], xlabel='RLD $(cm/cm^3)$')
plt.setp(axes[:, 0], ylabel='Depth $(cm)$')
plt.legend(np.asarray(legend_lst),loc="lower center", bbox_to_anchor=(-0.8, -0.5), ncol= 8)
fig.subplots_adjust()
plt.savefig("results/rld_plot.png",dpi=300, bbox_inches='tight')
plt.show() #|\label{l2_3:plotend}|