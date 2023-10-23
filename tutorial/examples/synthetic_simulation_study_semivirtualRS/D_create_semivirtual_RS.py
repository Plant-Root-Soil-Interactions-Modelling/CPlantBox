import sys
sys.path.append("../../../")
sys.path.append("../../../src/");
sys.path.append("../../../gui/viewer/");
import plantbox as pb
from functional.xylem_flux import XylemFluxPython 
import visualisation.vtk_plot as vp
import visualisation.vtk_tools as vt
from viewer_data import ViewerDataModel
from visualisation.vtk_tools import *
import rsml.rsml_writer
import rsml.rsml_reader as rsml
import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import os
import math
import random
#########################################################
def optimize_rootnum(rn0):
#root number is optimized

    ln = l_2lat / (rn0+rn3-len(idx))
    lpred = 0

    for i in range(0,len(polylines)):

        if int(funcs["subType"][i][0]) == 2:
            leng_cum = 0 # branch length up to the next strtnd 
            leng_tot = 0 #total branch length

            #compute total polyline length 
            d = np.diff(polylines[i], axis=0)
            polylen = np.sum(np.sqrt((d ** 2).sum(axis=1)))

            #print(len(polylines[i]))
            dummy = 0
            for j in range(0,len(polylines[i])-1):
                a = np.around(np.asarray(polylines[i][j]), decimals = 2) 
                b = np.around(np.asarray(polylines[i][j+1]), decimals = 2)

                if any((strtnd3[:]==a).all(1)): #if there is already a 3rd order lat
                    leng_cum = 0
                leng_tot = leng_tot + np.linalg.norm(a-b)
                if leng_tot>lb:
                    if dummy == 0:
                        lpred = lpred + (RSage+1-funcs["creationTime"][i][j])*r_virt
                    dummy = dummy +1
                    leng_cum = leng_cum + np.linalg.norm(a-b)

                #include angles & age & connections
                if leng_cum>ln: 
                    #find envisaged length of virtual segment
                    lpred = lpred + (RSage+1-funcs["creationTime"][i][j])*r_virt
                    leng_cum = 0

                if leng_tot>(polylen):
                    break

    err = (abs(misslen-lpred)/misslen) #relative error
    #print('There is still ' +  str(misslen-lpred) + ' cm root length missing!')
    return err


#########################################################

path  = 'results/'

#defined values 
fname = 'Faba_day10' #to be adjusted 
misslen = 97.47
rad = 0.025 #to be adjusted 
strtp = 10 #number of starting points for optimization 
lb = 1. #to be adjusted
la = 2. #to be adjusted
lmax = 4 #to be adjusted
theta = 1.423316 #to be adjusted
r_virt = 0.823 #to be adjusted
width = 2.8
depth = 19
bigcyl = pb.SDF_PlantContainer(width-0.05, width-0.05, depth-0.05, False)
RSage = 10

#read in the rsml scaffold 
rs = XylemFluxPython.read_rsml(path+ 'REC_' + fname +".rsml")
polylines, props, funcs, _ = rsml.read_rsml(path + 'REC_' + fname +".rsml")

#query nodes, segs 
ana = pb.SegmentAnalyser(rs)
nodes = np_convert(ana.nodes)
nodes = np.around(nodes, decimals=2)
segs = np_convert(ana.segments)

#store all original measures 
nodes_orig = nodes
segs_orig = segs
radius_orig = np.array(ana.getParameter("radius"))
length_orig = np.array(ana.getParameter("length"))
cts_orig = np.array(ana.getParameter("creationTime"))
type_orig = np.array(ana.getParameter("subType"), dtype = int)

#find out how many 3rd order lats are present and where they are located
typ = []
for sublist in funcs["subType"]:
    typ.append(sublist[0])
typ = np.asarray(typ)
rn1 = np.sum(typ==1)
rn2 = np.sum(typ==2)
rn3 = np.sum(typ==3)

strtnd3_ = np.zeros((rn3,3))
k = 0
for i in range(0,len(polylines)):
    if int(funcs["subType"][i][0]) == 3:
        strtnd3_[k,:] = np.around(polylines[i][0], decimals = 2) 
        k = k+1

#find the node before strtnd3_ (real start node) 
strtnd3 = np.zeros((len(strtnd3_),3))
for i in range(0, len(strtnd3_)): 
    idx3 = int(np.argwhere(np.all(nodes==strtnd3_[i,:],axis = 1)))
    idx3_seg = (np.where(segs[:,1]==idx3))[0]
    strtnd3[i,:] = nodes[segs[int(idx3_seg),0]]

#compute total branch lengths
brlength = []
for i in range(0,len(polylines)):
    d = np.diff(polylines[i], axis=0)
    brlength.append(np.sum(np.sqrt((d ** 2).sum(axis=1))))
brlength = np.asarray(brlength)
idx = np.where((typ == 2) & (brlength >=(la+lb)))[0]
l_2lat = np.sum(brlength[idx]-(lb)) #branching zone

#root number is optimized
rn0_ = np.zeros((strtp))
for i in range(1,len(rn0_)+1): 
    rn0_[i-1] = misslen / (lmax/i)

if l_2lat>0:
    best_fun = float('inf')
    er = np.zeros((len(rn0_)))
    rn_ = np.zeros((len(rn0_)))
    #bnds = np.hstack([0, rn0_])
    for n, rn0 in enumerate(rn0_):
        #print("initial guess root number = ", rn0)
        res_  = minimize(optimize_rootnum, rn0, method='Nelder-Mead') #optimize root number
        er[n] = res_.fun
        rn_[n] = res_.x[0]
    best = np.argmin(er)
    rn = rn_[best]
    ln = l_2lat / (rn+rn3-len(idx))
    print("optimized root number  = ", np.around(rn))
    #print(res.fun) 
else:
    ln = 10000 #to avoid any branches 
print('The inter-branch distance of the virtual roots is ', ln, ' cm')

#compute points of origin of missing roots  
strtnd = []; strtconn = []; strtang = []; strtage = []; lpred = 0; leng_tot_new = 0

for i in range(0,len(polylines)):

    if int(funcs["subType"][i][0]) == 2:
        leng_cum = 0 # branch length up to the next strtnd 
        leng_tot = 0 #total branch length 

        #compute total polyline length 
        d = np.diff(polylines[i], axis=0)
        polylen = np.sum(np.sqrt((d ** 2).sum(axis=1)))

        dummy = 0
        for j in range(0,len(polylines[i])-1):

            #check is max length is already reached and leave if yes 
            if leng_tot_new>misslen:
                print('missing root length is reached!')
                break

            a = np.around(np.asarray(polylines[i][j]), decimals = 2) 
            b = np.around(np.asarray(polylines[i][j+1]), decimals = 2)

            leng_tot = leng_tot + np.linalg.norm(a-b)
            if any((strtnd3[:]==a).all(1)): #check if there is already a 3rd order lat
                leng_cum = 0
                print(a, 'there is already a 3rd order lat!')

            #for angle computation
            c = np.array([1,1,0])/np.linalg.norm([1,1,0])
            d = (b-a)/np.linalg.norm(a-b)
            plunge = np.arccos(np.clip(np.dot(c,d), -1.0, 1.0))

            if leng_tot>=(lb):
                if dummy == 0:
                    strtnd.append(polylines[i][j])
                    strtang.append(plunge)
                    strtage.append(funcs["creationTime"][i][j])
                    #find envisaged length of virtual segment
                    lpred = lpred + (RSage+1-funcs["creationTime"][i][j])*r_virt

                    #find parent seg from parent node
                    polylines[i][j] = np.around(polylines[i][j], decimals=2)
                    testconn = np.argwhere(np.all(nodes==polylines[i][j],axis = 1))
                    strtconn.append(testconn[0][0])
                    leng_cum = 0

                    leng_tot_new = leng_tot_new + (RSage+1-funcs["creationTime"][i][j])*r_virt
                dummy = dummy+1
                leng_cum = leng_cum + np.linalg.norm(a-b)

            #include angles & age & connections
            if leng_cum>=ln: 
                strtnd.append(polylines[i][j])
                strtang.append(plunge)
                strtage.append(funcs["creationTime"][i][j])
                #find envisaged length of virtual segment
                lpred = lpred + (RSage+1-funcs["creationTime"][i][j])*r_virt

                #find parent seg from parent node
                polylines[i][j] = np.around(polylines[i][j], decimals=2)
                testconn = np.argwhere(np.all(nodes==polylines[i][j],axis = 1))
                strtconn.append(testconn[0][0])
                leng_cum = 0

                leng_tot_new = leng_tot_new + (RSage+1-funcs["creationTime"][i][j])*r_virt
            #print('total added length'+str(leng_tot_new) +'and total added rn = '+ str(len(strtnd)))
            if leng_tot>=(polylen):
                break

print('Total virtually added root length: ' + str(leng_tot_new))
print('Missing root length: ' + str(misslen))
            
#get RS params
path_rs = "../../../modelparameter/structural/rootsystem/"
name = "virtual"

# Initialize root systems
allRS = []
for i in range(0, len(strtnd)):
    rs = pb.RootSystem()
    rs.readParameters(path_rs + name + ".xml") 
    rs.getRootSystemParameter().seedPos = pb.Vector3d(strtnd[i][0], strtnd[i][1], strtnd[i][2])  # cm
    p = rs.getRootRandomParameter(1)
    p.theta = strtang[i]+p.theta #set theta dependent on parent
    p.lmax = lmax 
    p.a = rad
    rs.setSeed(i*random.randint(1,100))
    rs.setGeometry(bigcyl)
    rs.initialize()
    allRS.append(rs)


# Simulate
i = 0
typeVir = []; radiiVir = []; ctsVir = []; lengthVir = []; nodesVir = np.empty((0,3), int); segsVir = np.empty((0,2), int)
maxconn = len(nodes)-1

Leng = []
for rs in allRS:
   rs.simulate(RSage+1-strtage[i], True)
   r = np.asarray(rs.getParameter("radius"))
   t = np.ones((len(r)))*3        #set type to 3
   ct = np.asarray(rs.getParameter("creationTime"))
   ct = ct+strtage[i]                 #include start age of segment
   leng = np.asarray(rs.getParameter("length"))
   nodes = np.array(list(map(np.array,rs.getNodes())))
   seg = np.array(list(map(np.array,rs.getSegments())))

   lengi = np.zeros(len(seg))
   for q in range(0,len(seg)):
       a = nodes[seg[q,0],:]
       b = nodes[seg[q,1],:]
       lengi[q] = np.linalg.norm(a-b)
   Leng.append(sum(lengi))

   
   if len(seg):
       seg = seg+maxconn
       a = len(seg)+1
       seg_ = np.array([strtconn[i], seg[0,0]])  #include conencting segment
       seg = np.reshape(np.append(seg_,seg), (a,2))  
       maxconn = np.amax(seg)       #get new segment numbers for the rest of the segs


       nodesVir=np.append(nodesVir,nodes[1:,:],axis=0)
       segsVir= np.append(segsVir,seg,axis = 0)
       radiiVir = np.append(radiiVir, np.ones(len(seg))*r, axis=0)
       typeVir= np.append(typeVir, np.ones(len(seg))*t, axis = 0)
       ctsVir = np.append(ctsVir, np.ones(len(seg))*ct, axis = 0)
       lengthVir = np.append(lengthVir, np.ones(len(seg))*leng, axis = 0)
   i+=1


##connect original and virtual roots
radius = np.concatenate((radius_orig,radiiVir))
subtype = np.concatenate((type_orig, typeVir))
cts = np.concatenate((cts_orig,ctsVir))
length = np.concatenate((length_orig,lengthVir))
nodes = np.concatenate((nodes_orig,nodesVir),axis = 0)
segs = np.concatenate((segs_orig, segsVir), axis = 0)

radius = radius.tolist()
nodes = nodes.tolist()
cts = cts.tolist()
segs = segs.tolist()
subtype = subtype.tolist()
length = length.tolist()

nodes = [pb.Vector3d(n[0], n[1], n[2]) for n in nodes]
segs = [pb.Vector2i(s[0], s[1]) for s in segs]

ana = pb.SegmentAnalyser(nodes, segs, cts, radius)
ana.addData("radius", radius)
ana.addData("creationTime", cts)
ana.addData("subType", subtype)
ana.addData("length", length)
ana.write(path + "semivirt_"+fname+".vtp", ["radius", "subType", "creationTime", "length"] ) 

#write rsml
pd = vt.read_vtp(path + "semivirt_"+fname +'.vtp')
vt.write_rsml(path + "semivirt_"+fname + '.rsml', pd,0,[],[1])
print('done')
