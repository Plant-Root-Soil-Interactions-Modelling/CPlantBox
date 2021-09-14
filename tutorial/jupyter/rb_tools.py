import os
import sys
import numpy as np
sys.path.append("../..")
import plantbox as pb #CPlantBox Python Binding

#
# Auxiliary functions that could be moved to py_rootbox
#
def v2a(vd): # pb.std_vector_double_ to numpy array    
    l = np.zeros((len(vd),1)) 
    for i in range(0,len(vd)):
        l[i] = vd[i]
    return l

def v2ai(vd): # pb.std_vector_int_ to numpy int array    
    l = np.zeros(len(vd),dtype=np.int) 
    for i in range(0,len(vd)):
        l[i] = vd[i]
    return l

def a2v(a): #  numpy array to pb.std_vector_double
    l = pb.std_vector_double_()
    for d in a:
        l.append(d)
    return l

def a2i(a): #  numpy array to pb.std_vector_int
    l = pb.std_vector_int_()
    for i in a:
        l.append(i)
    return l
    
def vv2a(vd): # pb.std_vector_Vector3_ to numpy array
    N  = len(vd)
    l = np.zeros((N,3)) 
    for i in range(0,N):
        l[i,:] = [vd[i].x,vd[i].y,vd[i].z]
    return l

def seg2a(seg): # pb.std_vector_Vector2i_ to numpy array
    Ns = len(seg)
    seg_ = np.zeros((Ns,2),dtype = np.uint32)
    for i in range(0,Ns):
        seg_[i,:] = np.array([seg[i].x, seg[i].y])
    return seg_

def nodes2seg(nodes,seg,data): # node data to segment data 
    Ns = seg.shape[0]
    data_ = np.zeros(Ns)
    for i in range(0,Ns):
        n1 = seg[i,0]
        n2 = seg[i,1]
        data_[i] = 0.5*(data[n1]+data[n2])
    return data_
                
def z2i(z,n): # maps z to equidistant mesh
    i = int(round((abs(z)/100)*n))  
    return min(max(i,0),n-1) 

def plotRSinit(ax):
    ax.clear()
    ax.set_xlim(-0.2, 0.2)
    ax.set_ylim(-0.2, 0.2)
    ax.set_zlim(-0.35, 0)
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.set_zlabel("z [m]")         

def plotRS(ax,seg,nodes): # plots the root system as line plot (too slow)
    plotRSinit(ax)
    scale = 0.01;
    for s in seg:
        n1 = v2v(nodes[s.x])*scale
        n2 = v2v(nodes[s.y])*scale
        ax.plot((n1[0],n2[0]),(n1[1],n2[1]),(n1[2],n2[2]),'k-')
        
def plotRSscatter(ax,nodes): # plots the root system nodes (rather slow)
    plotRSinit(ax)
    scale = 0.01
    n = vv2a(nodes)
    ax.scatter(n[:,0]*scale,n[:,1]*scale,n[:,2]*scale)
    ax.set_title("Root tips")
        
def parameterPath(): # works only if everything is located in folder CPlantBox
    cwd = os.getcwd()
    i = cwd.index("CPlantBox"+os.sep)
    return cwd[0:i] + "CPlantBox"+os.sep+"modelparameter"+os.sep            
    
    