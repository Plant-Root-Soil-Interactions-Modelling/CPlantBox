"""
growing plant example, use of AnimateRoots 

TODO we could add how to create avi form png (e.g. on linux), and remove files again
"""
import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import visualisation.vtk_plot as vp

import numpy as np
from structural.Plant import PlantPython
import matplotlib.pyplot as plt

def get3Dshape(plant,title_ = 'wine', saveOnly = True):        
    orgs_ = plant.getOrgans(2)
    orgs_a = [1 for org in orgs_ if org.isAlive() ]
    orgs_al = [1 for org in orgs_ if (org.isAlive() and not org.getParameter('is_fine_root')) ]
    ana = pb.SegmentAnalyser(plant) 
    segOs = plant.getSegmentOrigins()
    
    #vp.plot_roots(ana, "id")

    '''
    Lignification status
    '''
    lignification = [segO.lignificationStatus() for segO in segOs]
    '''
    Survival
    '''
    aliveSegs = [segO.isAlive() for segO in segOs]
    
    print('alive nodes:', sum(aliveSegs), 'alive organs',sum(orgs_a), 'alive long lived roots',sum(orgs_al))
    ana.addData('alive', aliveSegs)
    ana.addData('lignification', lignification)
    ana.filter('alive', 1)
    vp.plot_roots(ana, "subType",p_names = ['lignification',"creationTime","id"] , win_title = title_, render = not saveOnly)
    
doAnim = False

plant = PlantPython(1)

# Open plant and root parameter from a file
path = "../../modelparameter/structural/rootsystem/"
name = "wine_Fichtl"

plant.readParameters(path + name + ".xml")


ps = plant.getOrganRandomParameter(pb.seed)[0]
ps.Lmax_unsuberized = 5.
ps.Lmax_suberized = 10.
ps.delayDefinition = 4




for ii, pp in enumerate(plant.getOrganRandomParameter(2)):
    if ii < 2:
        pp.ldelay  = 0*1e6
        pp.ldelays = 1e6*50#*5
    else:
        pp.ldelay  = 0*1e6
        pp.ldelays = 1e6*50#/2.
        
allRRP = plant.getOrganRandomParameter(2)
allRRP[2].a_gr = 0.085/1e6
allRRP[3].a_gr = 0.085/1e6
allRRP[3].k_survive = 2.17
allRRP[3].lambda_survive = 4.55

for newpId in range(3,6):
    rrpA = allRRP[newpId].copy(plant)
    
    rrpA.subType += 1
    rrpA.successorST =  np.array(rrpA.successorST) + 1
    rrpA.successorP[0][0] *= 0.5
    rrpA.successorP[0][1] *= 0.75
    #rrpA.successorP =  np.array(rrpA.successorP) *0.5
    rrpA.a_gr = 0.085/1e6
    #rrpB.successorP =  np.array(rrpB.successorP) /10.
    allRRP.append(rrpA)
    
allRRP[6].successorST = [[11]]
allRRP[6].successorP = [[allRRP[6].successorP[0][1]]] # [[pp1],[pp2]]

pplats = np.array([[0.35]])#,[0.35]])
for newpId in range(6,11):
    rrpA = allRRP[newpId].copy(plant)
    
    rrpA.subType += 1
    snext = rrpA.subType + 1
    rrpA.successorST =  [[snext]] # [[snext],[snext]]
    #pp = rrpA.successorP[0][0] * 0.5
    rrpA.successorP =  pplats
    #rrpB.successorP =  np.array(rrpB.successorP) /10.
    allRRP.append(rrpA)
    rrpA.is_fine_root = True
    rrpA.lmax = ps.Lmax_suberized
    
    pplats *= 0.5

allRRP[-1].successorST = []

for newp in allRRP:
    plant.setOrganRandomParameter(newp)
   

Main_beta = 3.28 * 10e-5
Sub_beta = 2.16 * 10e-5
SubSub_beta = 2.00 * 10e-5

p2 = plant.getOrganRandomParameter(2, 2)
p2.r = p2.lmax * Main_beta

p2.tropismT = 7  # mix of planar and gravitropism
p2.tropismS = 0.2
p2.tropismN = 0.5  
p2.tropismW1 = 0.4
p2.tropismW2 = 0.6


p3 = plant.getOrganRandomParameter(2, 3)
p3.r = p3.lmax * Sub_beta

for subindx in [4,5,6]:
    p4 = plant.getOrganRandomParameter(2, subindx)
    p4.r = p4.lmax * SubSub_beta
    
    
# for ii, pp in enumerate(plant.getOrganRandomParameter(2)):
#    print(pp)
    
'''
check Tropisms
'''
# for pp in plant.getOrganRandomParameter(2):
    # pp.tropismT = 7  # mix of planar and gravitropism
    # pp.tropismS = 1
    # pp.tropismN = 10
    # pp.tropismW1 = 1.
    # pp.tropismW2 = 0.


'''
check inputs
'''
#for ii, pp in enumerate(plant.getOrganRandomParameter(2)):
#    print(pp)
    
plant.initialize_static(path + "B-23_Fichtl.rsml", [0, 1])  # 0 is shoot, 1 are static roots

# the static laterals 2 and 3 are replaced with the growing lateral 2
ld = plant.set_identical_laterals([0, 1], [2, 3], 2)
plant.initialize_static_laterals()

plt.hist(np.array(ld)/1e6, density = False, bins=30)
plt.title("Creation time of the main roots")
plt.show()

dt = 1e6 # ~1 yr
N = 50
min_ = np.array([-20, -20, -50])
max_ = np.array([20, 20, 30.])



if doAnim:
    anim = vp.AnimateRoots(plant)
    anim.min = min_
    anim.max = max_
    anim.res = [1, 1, 1]
    anim.file = "results/example_plant"
    anim.avi_name = "results/example_"
    anim.plant = True
    anim.start()
    
for i in range(0, N):
    if i > 0:
        plant.survivalTest()
    plant.simulate(dt, False)
    print('time', dt * i /1e6)
    get3Dshape(plant,title_ = 'wine'+str(i), saveOnly = True)
        
    if doAnim:
        anim.root_name = "organType"
        anim.update()
        

orgs_ = plant.getOrgans(2)
orgs_long = [org for org in orgs_ if org.param().subType in range(2,7)]
orgs_fine = [org for org in orgs_ if org.param().subType in range(7,12)]

'''
Age distribution at last time step
'''
ages_long = np.array([org.getAge()/1e6 for org in orgs_long])
subType_long = np.array([org.param().subType for org in orgs_long])
print(set(subType_long))
ages_long = [ages_long[np.isin(subType_long, hh)] for hh in [[2],[3],[4,5,6]]]

plt.hist(ages_long, 
         bins=30, 
         density=False, 
         stacked = True,
         #histtype='stepfilled',
         #alpha=0.5,          
         edgecolor='black',
            label=['main','sub', 'subsub'] )
plt.title("Living and dead roots age distribution at last time step")
plt.legend()
plt.show()

'''
secondary growth
'''
for org in orgs_long:
    age_segs = (org.getAge() + org.getNodeCT(0) - np.array([org.getNodeCT(ii) for ii in range(1,org.getNumberOfNodes())]))/1e6
    plt.plot(age_segs,org.getRadii())
plt.title("Radius vs time for long lived roots")
plt.show()


'''
Survival rate
'''
rtl_winters = np.array([org.getParameter('rlt_winter')/1e6 for org in orgs_long])
rtl_winters = [rtl_winters[np.isin(subType_long, hh)] for hh in [[2],[3],[4,5,6]]]
#maxT = int(max(rtl_winters)) + 1
outout = [[],[],[]]
for idrt, RootType in enumerate(['main','sub', 'subsub']):
    if len(rtl_winters[idrt] > 0):
        maxT = int(max(rtl_winters[idrt])) + 1
        for yr in range(maxT):
            ratioAlive = sum(rtl_winters[idrt] > yr)/len(rtl_winters[idrt])
            outout[idrt].append(ratioAlive)
        plt.scatter([yr for yr in range(maxT)],outout[idrt], label=RootType)
plt.title("Survival probability")
plt.legend()
plt.show()


'''
Creation time
'''
creationTime_long = np.array([org.getParameter('creationTime')/1e6 for org in orgs_long])
creationTime_short = np.array([org.getParameter('creationTime')/1e6 for org in orgs_fine])

plt.figure()
plt.hist([creationTime_long,creationTime_short], 
            stacked=True, density = False, bins=30, edgecolor='black',
            label=[ 'long-lived roots','fine roots'] )
plt.title("Root emergence time")
plt.legend()
plt.show()
