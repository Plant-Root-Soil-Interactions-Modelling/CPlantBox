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

'''
length vs time init
'''
# def getlen1(BEDD_tot,beta,L_max):
    # return L_max *(1-np.exp(-beta * BEDD_tot))

# def getlen2(BEDD_tot,r,L_max):
    # return L_max *(1-np.exp(-r/L_max * BEDD_tot))
    
# Mainbeta = 3.28 * 10e-5
# Subbeta = 2.16 * 10e-5
# SubSubbeta  = 2.00 * 10e-5
# L_max = 96
# r = L_max * Mainbeta

# testmain = [getlen1(btot,Mainbeta,L_max) for btot in range(int(1225))]
# testmain2 = [getlen2(btot,r,L_max) for btot in range(int(1225))]
# btots = [btot for btot in range(int(1225))]
# plt.plot(btots, testmain)
# #plt.show()
# plt.plot(btots, testmain2)
# plt.show()

def get3Dshape(plant,title_ = 'wine', saveOnly = True):        
    orgs_ = plant.getOrgans(2)
    orgs_a = [1 for org in orgs_ if org.isAlive() ]
    orgs_al = [1 for org in orgs_ if (org.isAlive() and not org.getParameter('is_fine_root')) ]
    ana = pb.SegmentAnalyser(plant) 
    segOs = plant.getSegmentOrigins()
    
    #vp.plot_roots(ana, "id")

    '''
    Lignification status, Survival, fine roots
    '''
    lignification = [segO.lignificationStatus() for segO in segOs]
    aliveSegs = [segO.isAlive() for segO in segOs]
    is_fine_root = [segO.getParameter('is_fine_root') for segO in segOs]
    
    print('alive nodes:', sum(aliveSegs), 'alive organs',sum(orgs_a), 'alive long lived roots',sum(orgs_al))
    ana.addData('alive', aliveSegs)
    ana.addData('lignification', lignification)
    ana.addData('is_fine_root', is_fine_root)
    ana.filter('alive', 1)
    vp.plot_roots(ana, "subType",p_names = ['lignification','is_fine_root',"creationTime","id"] , win_title = title_, render = not saveOnly)
    
doAnim = False

soilSpace = pb.SDF_PlantContainer(np.inf, np.inf,  np.inf, True)  # to avoid root growing aboveground
plant = PlantPython(1)

# Open plant and root parameter from a file
path = "../../modelparameter/structural/rootsystem/"
name = "wine_Fichtl"

plant.readParameters(path + name + ".xml")

plant.setGeometry(soilSpace) 

ps = plant.getOrganRandomParameter(pb.seed)[0]
ps.Lmax_unsuberized = 5.
ps.Lmax_suberized = 10.
ps.delayDefinition = 4

yr_to_BEDD = 1225


for ii, pp in enumerate(plant.getOrganRandomParameter(2)):
    if ii == 0:
        pp.ldelay  = 0*yr_to_BEDD
        pp.ldelays = yr_to_BEDD*200#*5
    elif ii < 2:
        pp.ldelay  = 0*yr_to_BEDD
        pp.ldelays = yr_to_BEDD*40#*5
    else:
        pp.ldelay  = 0*yr_to_BEDD
        pp.ldelays = yr_to_BEDD*40#/2.
    pp.a_gr =  0.083/2/yr_to_BEDD

        
allRRP = plant.getOrganRandomParameter(2)


for newpId in range(3,6):
    rrpA = allRRP[newpId].copy(plant)
    
    rrpA.subType += 1
    rrpA.successorST =  np.array(rrpA.successorST) + 1
    #rrpA.successorP[0][0] *= 0.5
    #rrpA.successorP[0][1] *= 0.75
    #rrpA.successorP =  np.array(rrpA.successorP) *0.5
    rrpA.a_gr = 0.083/2/yr_to_BEDD
    #rrpB.successorP =  np.array(rrpB.successorP) /10.
    allRRP.append(rrpA)
    
allRRP[4].k_survive = 2.17
allRRP[4].lambda_survive = 4.55

allRRP[6].successorST = [[11]]
allRRP[6].successorP = [[allRRP[6].successorP[0][1]]] # [[pp1],[pp2]]

pplats = np.array([[0.5]])#,[0.35]])
for newpId in range(6,11):
    rrpA = allRRP[newpId].copy(plant)
    
    rrpA.subType += 1
    snext = rrpA.subType + 1
    rrpA.successorST =  [[snext]] # [[snext],[snext]]
    #pp = rrpA.successorP[0][0] * 0.5
    rrpA.successorP =  pplats
    #rrpB.successorP =  np.array(rrpB.successorP) /10.
    rrpA.is_fine_root = True
    rrpA.lmax = 5. #ps.Lmax_suberized
    allRRP.append(rrpA)
    
    pplats *= 0.75

allRRP[-1].successorST = []

for newp in allRRP:
    #newp.tropismS = 0.8
    plant.setOrganRandomParameter(newp)
   

Main_beta = 3.28 * 10e-5
Sub_beta = 2.16 * 10e-5
SubSub_beta = 2.00 * 10e-5

p2 = plant.getOrganRandomParameter(2, 2)
p2.r = p2.lmax * Main_beta

p2.tropismT = 7  # mix of planar and gravitropism
p2.tropismN = 2  
p2.tropismW1 = 0.85 #0.4 # gravitropism
p2.tropismW2 = 0.15 #0.6 # plagiotropism


p3 = plant.getOrganRandomParameter(2, 3)
p3.r = p3.lmax * Sub_beta

for subindx in range(4,11):
    p4 = plant.getOrganRandomParameter(2, subindx)
    p4.r = p3.r # p4.lmax * SubSub_beta
    
    
#for ii, pp in enumerate(plant.getOrganRandomParameter(2)):
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
    
    
    
'''
start simulation
'''
plant.initialize_static(path + "B-23_Fichtl.rsml", [0, 1])  # 0 is shoot, 1 are static roots

# the static laterals 2 and 3 are replaced with the growing lateral 2
ld, ld1 = plant.set_identical_laterals([0, 1], [1, 2, 3], 2)
# plt.hist(np.array(ld1)/yr_to_BEDD, density = False, bins=30)
# plt.title("Creation time of the main roots")
# plt.show()
# plt.hist(np.array(ld)/yr_to_BEDD, density = False, bins=30)
# plt.title("Creation time of the main roots")
# plt.show()

plant.initialize_static_laterals()
plant.betaN = 5000
  

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
    
# N1 = 365
# dt = yr_to_BEDD/N1 # ~1 d

# all_lengths = []
# all_ages = []
# all_subtypes = []

# for i in range(0, N1):
    # plant.simulate(dt, False)
    # print('time', dt * i /yr_to_BEDD, end=", ")
    # get3Dshape(plant,title_ = 'wine'+str(i), saveOnly = True)
        
    # orgs_ = plant.getOrgans(2)
    # orgs_ = [org for org in orgs_ if org.param().subType > 0] # ignore the static roots
    # all_ages.append(np.array([org.getAge() for org in orgs_]))
    # all_lengths.append(np.array([org.getLength() for org in orgs_]))
    # all_subtypes.append(np.array([org.param().subType for org in orgs_]))
    
    # if doAnim:
        # anim.root_name = "organType"
        # anim.update()
        
# '''
# Length vs age
# '''
# all_ages = np.array([item for sublist in all_ages for item in sublist])
# all_lengths = np.array([item for sublist in all_lengths for item in sublist])
# all_subtypes = np.array([item for sublist in all_subtypes for item in sublist])
# for stroots in set(all_subtypes):
    # age_ = all_ages[all_subtypes == stroots]
    # len_ = all_lengths[all_subtypes == stroots]
    # plt.scatter(age_,len_, label = stroots)
# plt.title("Length vs age")
# plt.legend()
# plt.show()


all_lengths = []
all_ages = []
all_subtypes = []
all_alive = []
dt = yr_to_BEDD  # ~1 yr
N = 49

# orrs = plant.getOrgans(-1)[0]
# print(orrs.organType(),orrs.param().subType,orrs.getAge(),orrs.getLength(),
            # orrs.getRootRandomParameter().a_gr, orrs.getRadius(0))
            
# orr = plant.getOrgans(2)[0]
# print(orr.param().subType,orr.getNumberOfChildren(), orr.getParameter('age')/yr_to_BEDD,orr.getParameter('rlt_winter')/yr_to_BEDD  )
# for orrid in range(orr.getNumberOfChildren()):
    # orrk = orr.getChild(orrid)
    # print(orrk.param().subType,orrk.getNumberOfChildren(), orrk.getParameter('age')/yr_to_BEDD,orrk.getParameter('rlt_winter')/yr_to_BEDD  )

# raise Exception
get3Dshape(plant,title_ = 'wine0', saveOnly = True)
for i in range(0, N):
    plant.survivalTest()
    plant.simulate(dt, False)
    print('time', dt * i /yr_to_BEDD, end=", ")
    get3Dshape(plant,title_ = 'wine'+str(i+1), saveOnly = True) #+N1
        
    orgs_ = plant.getOrgans(2)
    #orgs_ = [org for org in orgs_ if org.param().subType > 0] # ignore the static roots
    all_ages.append(np.array([org.getAge()/yr_to_BEDD for org in orgs_]))
    all_lengths.append(np.array([org.getLength() for org in orgs_]))
    all_alive.append(np.array([org.isAlive() for org in orgs_]))
    all_subtypes.append(np.array([org.param().subType for org in orgs_]))
    
    zzz = np.array([xyz[2] for xyz in plant.get_nodes()])
    print('zzz[zzz > 0]',zzz[zzz > 0])
        
    # orrs = plant.getOrgans(-1)[0]
    # print(orrs.getId(),orrs.organType(),orrs.param().subType,orrs.getAge(),orrs.getLength(),
                # orrs.getRootRandomParameter().a_gr, orrs.getRadius(0),
                # # orrs.getNodeCT(0))
    # if i == 2:
        # raise Exception
    
    if doAnim:
        anim.root_name = "subType"
        anim.update()
    

#orgs_ = plant.getOrgans(2)
orgs_long = [org for org in orgs_ if org.param().subType in range(2,7)]
orgs_fine = [org for org in orgs_ if org.param().subType in range(7,12)]

'''
Presence (creation time vs death)
'''
#all_alive = np.array([item for sublist in all_alive for item in sublist])
#all_subtypes_ = np.array([item for sublist in all_subtypes for item in sublist])
alive_long = [sum(all_alive[i][np.isin(all_subtypes[i],[2,3,4,5,6])]) for i in range(len(all_alive))]
alive_short = [sum(all_alive[i][np.isin(all_subtypes[i],[7,8,9,10,11])]) for i in range(len(all_alive))]
#alive_short = all_alive[np.isin(all_subtypes_,[2,3,4,5,6])]
#alive_long = all_alive[np.isin(all_subtypes_,[7,8,9,10,11])]
year_ = [i for i in range(len(all_alive))]
plt.figure()
plt.scatter(year_,alive_long, label= 'long-lived roots')
plt.scatter(year_,alive_short, label='fine roots')
#plt.hist([alive_long,alive_short], 
#            stacked=True, density = False, bins=30, edgecolor='black',
#            label=[ 'long-lived roots','fine roots'] )
plt.title("Root count")
plt.legend()
plt.show()

'''
Length vs age
'''
all_ages_ = np.array([item for sublist in all_ages for item in sublist])
all_lengths_ = np.array([item for sublist in all_lengths for item in sublist])
all_subtypes_ = np.array([item for sublist in all_subtypes for item in sublist])
for stroots in set(all_subtypes_):
    age_ = all_ages_[all_subtypes_ == stroots]
    len_ = all_lengths_[all_subtypes_ == stroots]
    plt.scatter(age_,len_, label = stroots)
plt.title("Length vs age")
plt.legend()
plt.show()

'''
Age distribution at last time step
'''
ages_long_ = np.array([org.getAge()/yr_to_BEDD for org in orgs_long])
subType_long_ = np.array([org.param().subType for org in orgs_long])
subType_fine_ = np.array([org.param().subType for org in orgs_fine])
print('root subtypes obtained',set(subType_long_),set(subType_fine_))
ages_long_ = [ages_long_[np.isin(subType_long_, hh)] for hh in [[2],[3],[4,5,6]]]

plt.hist(ages_long_, 
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
    age_segs = (org.getAge() + org.getNodeCT(0) - np.array([org.getNodeCT(ii) for ii in range(1,org.getNumberOfNodes())]))/yr_to_BEDD
    plt.plot(age_segs,org.getRadii())
plt.title("Radius vs time for long lived roots")
plt.show()


'''
Survival rate
'''
rtl_winters = np.array([org.getParameter('rlt_winter')/yr_to_BEDD for org in orgs_long])
rtl_winters = [rtl_winters[np.isin(subType_long_, hh)] for hh in [[2],[3],[4,5,6]]]
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

import pickle
with open('./SurvivalRate.pkl','wb') as f:
    pickle.dump(rtl_winters,f, protocol=pickle.HIGHEST_PROTOCOL)

'''
Creation time
'''
creationTime_long = np.array([org.getParameter('creationTime')/yr_to_BEDD for org in orgs_long])
creationTime_short = np.array([org.getParameter('creationTime')/yr_to_BEDD for org in orgs_fine])

plt.figure()
plt.hist([creationTime_long,creationTime_short], 
            stacked=True, density = False, bins=30, edgecolor='black',
            label=[ 'long-lived roots','fine roots'] )
plt.title("Root emergence time")
plt.legend()
plt.show()

