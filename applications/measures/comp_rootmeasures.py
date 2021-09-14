import sys
sys.path.append("../..")
sys.path.append("../../src/python_modules/")
import time
import numpy as np
import plantbox as rb
import matplotlib
import matplotlib.pyplot as plt
from scipy.spatial import Delaunay
import math

def tetrahedron_volume(a, b, c, d):
    return np.abs(np.einsum('ij,ij->i', a-d, np.cross(b-d, c-d))) / 6

def convHull(pts):
    dt = Delaunay(pts)
    tets = dt.points[dt.simplices]
    vol = np.sum(tetrahedron_volume(tets[:, 0], tets[:, 1],tets[:, 2], tets[:, 3]))
    return vol

def uniqueIndexes(l):
    seen = set()
    res = []
    for i, n in enumerate(l):
        if n not in seen:
            res.append(i)
            seen.add(n)
    return res

#############################################################

path = 'dgf_files'

##
#read in dgf
##
name = 'DAP8.dgf'
with open(path+'/'+name) as f:
    content = f.readlines()

i = 3
coordlist = []
while ('#' in content[i])==False: 
    line = content[i]
    nums = line.split()
    nums = map(float,nums)
    nums = list(nums)
    coordlist.extend((nums))
    i += 1        
A = np.asarray(coordlist)
nodes = np.reshape(A, (int(len(A)/3), 3))

i = i+3
paramlist = []
while ('#' in content[i])==False: 
    line = content[i]
    nums = line.split()
    nums = map(float,nums)
    nums = list(nums)
    paramlist.extend((nums))
    i += 1        
B = np.asarray(paramlist)
params = np.reshape(B, (int(len(B)/10), 10))

##
#parameters...
##
x = nodes[2:-1,0]*100;
y = nodes[2:-1,1]*100;
z = nodes[2:-1,2]*100; 
leng = params[:,5]; 
rad = params[:,6]; 
etime = params[:,9]; 
branch = params[:,3]; 
prev = params[:,0]; 
order = params[:,2];

##
#compute max root age  
##
etime = np.ceil(etime) 
maxage = int(np.max(etime))


##
#preallocate root measures and compute them 
##
TRL = np.zeros((maxage,))
gr = np.zeros((maxage,))
convH = np.zeros((maxage,))
rootnum = np.zeros((maxage,))
num1lat = np.zeros((maxage,))
num2lat = np.zeros((maxage,))
RLD = np.zeros((maxage,))
HMD = np.zeros((maxage,))

dummy = 0
for i in range(1,maxage):
    idx_ = np.where(etime <= i)
    idx = idx_[0]
    TRL[i] = np.sum(leng[idx])
    gr[i] = TRL[i]-dummy
    convH[i] = convHull(nodes[idx,:]*100) #from m to cm
    rootnum[i] = len(np.unique(branch[idx]))
    idxbr = uniqueIndexes(branch[idx])
    num1lat[i] = np.sum(order[idxbr]==2)
    num2lat[i] = np.sum(order[idxbr]==3)
    RLD[i] = TRL[i]/convH[i]
    HMD[i] = (math.pi*RLD[i])**-0.5

plt.rcParams.update({'font.size': 16})
fig,axs = plt.subplots(2,3)
x = np.linspace(1,maxage,maxage)
axs[0,0].plot(x,TRL, 'r')
axs[0,0].set_xlabel('Time (days)')
axs[0,0].set_ylabel('Total root length (cm)')
axs[0,0].legend([name])

axs[0,1].plot(x,gr, 'r')
axs[0,1].set_xlabel('Time (days)')
axs[0,1].set_ylabel('Growth rate of \n the root system (cm $d^{-1})$')

axs[0,2].plot(x,convH, 'r')
axs[0,2].set_xlabel('Time (days)')
axs[0,2].set_ylabel('Convex hull $(cm^{3})$')

axs[1,0].plot(x,rootnum, 'r')
axs[1,0].set_xlabel('Time (days)')
axs[1,0].set_ylabel('Total number of root tips (-)')

axs[1,1].plot(x,RLD, 'r')
axs[1,1].set_ylim(0,3000)
axs[1,1].set_xlabel('Time (days)')
axs[1,1].set_ylabel('Root length \n  density (cm $cm^{-3})$')

axs[1,2].plot(x,HMD, 'r')
axs[1,2].set_xlabel('Time (days)')
axs[1,2].set_ylabel('Half mean distance (cm)')

plt.tight_layout()
plt.subplots_adjust(wspace=0.4,hspace=0.3)
plt.show()

#np.savez('data/rootmeasures.npz', x, TRL, gr, convH, rootnum, RLD, HMD, num1lat, num2lat) 

# root measures
# (1) total root length - TRL
# (2) growth rate of the root system (cm/d) - gr
# (3) convex hull -convH
# (4) total number of roots - rootnum
# (5) RLD (=total root length / convex hull) - RLD
# (6) HMD (=(pi*RLD)^-0.5) - HMD
# (7) number of 1st order laterals - num1lat
# (8) number of 2nd order laterals - num2lat
