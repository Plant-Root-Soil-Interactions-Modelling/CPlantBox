import sys
sys.path.append("../../")
sys.path.append("../../src")
sys.path.append("../../modelparameter/structural/rootsystem/")
import plantbox as pb
import os
import numpy as np
import plantbox as pb
import multiprocessing as mp
import time
import csv

#General parameters 
name = "Zea_mays_3_Postma_2011"
simtime = 120 #days
M = 8 # number of plants in one row
N = 6 # number of rows
distp = 12  # distance between the plants[cm] (P-P)
distr = 75 # distance between the rows[cm] (R-R)
l_field=M*distp # length of field 
w_field=N*distr # width of field 
depth = 130

#Geometry parameters for  rhizotubes 
tube_diam = 6.4 #tube diameter
tube_len = 200 #tube length 
y_ = [-25, -15, -5, 5, 15, 25] #y positions of rhizotubes 
x_ = [-25, -50, -75, -100, -125, -150] #x positions of rhizotubes 
incl_angle = 60 #inclination angle from the vertical in (°)

#Simulate rhizotubes
rhizotube_ = pb.SDF_PlantContainer(tube_diam/2, tube_diam/2, tube_len, False)  # a single rhizotube
rhizotube = pb.SDF_RotateTranslate(rhizotube_, incl_angle, pb.SDF_Axis.yaxis, pb.Vector3d(tube_len / 2, 0, 0))  # a single rhizotube with correct x position 
rhizotubes_ = []
for i in range(0, len(y_)):
    rhizotubes_.append(pb.SDF_RotateTranslate(rhizotube, pb.Vector3d(tube_len/2+x_[i], y_[i], 0)))
rhizotubes = pb.SDF_Union(rhizotubes_) 
rs = pb.Plant()
rs.setGeometry(rhizotubes)
if not os.path.exists('results_rhizotubes/'):
   os.makedirs('results_rhizotubes/')
rs.write("results_rhizotubes/rhizotubes.py")

# Simulate N*M root systems
allRS = []
for i in range(0, M):
    for j in range(0, N):
        rs = pb.Plant()
        rs.readParameters('../../modelparameter/structural/rootsystem/'+name + ".xml")
        params = rs.getOrganRandomParameter(pb.root)
        seed = rs.getOrganRandomParameter(pb.seed)[0]
        seed.seedPos = pb.Vector3d(distr * j-((N-1)*distr/2), distp * i-(distp*M/2), 0) 
        rs.initialize(False)
        rs.simulate(simtime, False)
        allRS.append(rs)
        if i+ j == 0:
            allAna = pb.SegmentAnalyser(rs)
        else:
            allAna.addSegments(rs) 
allAna.mapPeriodic(w_field, l_field) #periodic outer boundaries 
allAna.write("results_rhizotubes/RS_field.vtp")


