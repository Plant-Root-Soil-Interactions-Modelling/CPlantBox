import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb

import os
import pstats
from io import StringIO
from datetime import datetime, timedelta
import numpy as np
import random
from functional.plant_conductivities import init_conductivities
from functional.phloem_flux import PhloemFluxPython  # root system Python hybrid solver

######################
#
# plant
#
#######################
pl = pb.MappedPlant(seednum = 1)  # set seed
plFlow = PhloemFluxPython(pl,psiXylInit = -659.8,ciInit = 350e-6*0.5) 
path = "../../modelparameter/structural/plant/"
name = "Triticum_aestivum_test_2021"  # root growth model parameter file

pl.readParameters(path + name + ".xml")
simtime = 10
dt = simtime
steps = int(simtime/dt)
pl.initialize(verbose = False)
init_conductivities(plFlow)
for step in range(steps):
    print("\n\n\nstep n°", step)
    pl.simulate(dt, False)
    
soil_index = lambda x, y, z: int((z<=0)*2 -1)
pl.setSoilGrid(soil_index)


orgs = pl.getOrgans(2, False)# take all the roots
st_orgs = np.array([org.getParameter("subType") for org in orgs])
lenOrg = np.array([org.getLength(False) for org in orgs])   
globalSegId = np.array([ np.array(org.getNodeIds()[1:]) for org in orgs],dtype=object) -1#.reshape(-1)
a = np.array(pl.radii)# radius
st =np.array(pl.subTypes )# conductivities kr, kx 
ot = plFlow.get_organ_types()
dist = np.array(pl.distanceTip)

#print(globalSegId)
#print(plFlow.kr_f(-10, st, ot, seg_ind))
for org in [orgs[0]]:
    print('id', org.getId(), 'ot',org.organType(),'st',org.getParameter("subType"),'length',org.getParameter("length"))
    gSegId = np.array(org.getNodeIds()[1:])-1
    print('globalSegId',gSegId)
    print('dist',dist[gSegId])
    print('st',st[gSegId], 'ot',ot[gSegId])
    print([plFlow.kr_f_cpp(seg_ind, 0, st[seg_ind], ot[seg_ind]) for seg_ind in gSegId])
    print("\n")
