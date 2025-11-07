""""more complex geometries"""
import sys; sys.path.append("../../.."); sys.path.append("../../../src/")
import plantbox as pb
import visualisation.vtk_plot as vp
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cm as cm
import os
import matplotlib.colors

mpl.rcParams['mathtext.default'] = 'regular'
mpl.rcParams.update({'font.size': 20})
save = True

rs = pb.RootSystem()
path = "data/"
name = "virtual"
rs.readParameters(path + name + ".xml")

#circular Obstacles with 0.25cm, 0.35cm, 0.5cm, 1cm radius
rad = 0.5 #radius - can eb changed 
depth = 1 #depth of rhizotron box
S = 0.5 #horizontal distance between two obstacles 
box = pb.SDF_PlantBox(depth, 14*rad+5.5*S, 2*rad+7*((2*rad+S)**2-(rad+S/2)**2)**0.5)  # box
obstacle = pb.SDF_PlantContainer(rad, rad, depth, False)  # a single obstacle
rhizoX = pb.SDF_RotateTranslate(obstacle, 90, pb.SDF_Axis.yaxis, pb.Vector3d(depth / 2, 0, 0))

y1_= np.arange(rad, rad+5*(2*rad+S), (2*rad+S))-(14*rad+5.5*S)/2
y2_= np.arange((2*rad+0.5*S), (2*rad+0.5*S)+5*(2*rad+S), (2*rad+S))-(14*rad+5.5*S)/2
z_= -1*np.arange(rad,rad+7*((2*rad+S)**2-(rad+S/2)**2)**0.5,((2*rad+S)**2-(rad+S/2)**2)**0.5)

obstacles_ = []
for i in range(0, len(y1_)):
    for j in range(0, len(z_)):
        if j%2: 
            y_ = y1_
        else: 
            y_ = y2_
        v = pb.Vector3d(0, y_[i], z_[j])
        obstacles_.append(pb.SDF_RotateTranslate(rhizoX, v))

obstacles = pb.SDF_Union(obstacles_)
rs.setGeometry(obstacles)
if save: 
    rs.write("results/obstacles_rad"+str(int(rad*10))+".py")
obstacles_diff = pb.SDF_Difference(box, obstacles)
rs.setGeometry(obstacles_diff)

# Simulate
tropismN_ = np.linspace(1,4,4) 
tropismS_ = np.linspace(0.025, 0.2, 8) 
fig, ax = plt.subplots(len(tropismN_), len(tropismS_))
for i in range(0, len(tropismN_)): 
    for j in range(0, len(tropismS_)): 
        rs.getRootRandomParameter(1).tropismN=tropismN_[i]
        rs.getRootRandomParameter(1).tropismS=tropismS_[j]
        rs.getRootRandomParameter(1).dx=0.1
        rs.setSeed(1)
        rs.initialize()
        
        name = "Yao_exp_N"+str(tropismN_[i])+'_S'+str(np.around(tropismS_[j],3))
        if not os.path.exists('results/rad'+str(int(rad*10))+'/'+name):
            if save: 
                os.makedirs('results/rad'+str(int(rad*10))+'/'+name)
        
        for k in range(0,8): 
            rs.simulate(1)  # days
            if save: 
                rs.write('results/rad'+str(int(rad*10))+'/'+name+"/Yao_exp_N"+str(tropismN_[i])+'_S'+str(np.around(tropismS_[j],3))+'_day'+str(k)+".vtp")
        ana = pb.SegmentAnalyser(rs)
        segs = ana.segments
        nodes = ana.nodes
        for m, s in enumerate(segs):
            s1 = segs[m]
            n1, n2 = nodes[s1.x], nodes[s1.y]
            ax[i,j].plot([n1.y, n2.y], [n1.z, n2.z],color = 'r')
            
        ax[i,j].axis('equal')
        ax[i,j].set_xticks([])
        ax[i,j].set_yticks([])
        if i == (len(tropismN_)-1): 
            ax[i,j].set_xlabel('S = '+str(np.around(tropismS_[j],3)))
        if j == 0: 
            ax[i,j].set_ylabel('N = '+str(int(tropismN_[i])))

plt.show()