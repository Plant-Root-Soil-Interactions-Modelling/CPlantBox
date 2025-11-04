"""
Get the new parameters 
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
   
def run_sumlation(params):   
    doAnim = False

    soilSpace = pb.SDF_PlantContainer(np.inf, np.inf,  np.inf, True)  # to avoid root growing aboveground
    plant = PlantPython()

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
            pp.ldelays = yr_to_BEDD * params['ldelays0'] #200#*5
            pp.successorP[0][0] = params['successorP0'] 
            pp.successorNo[0] = params['successorNo0'] 
        else:
            pp.ldelay  = 0*yr_to_BEDD
            pp.ldelays = yr_to_BEDD * params['ldelays1']
            pp.successorP[0][0] = params['successorP10'] 
            pp.successorP[0][1] = params['successorP11'] 
            pp.successorNo[0] = params['successorNo0'] 
        pp.a_gr =  0.083/2/yr_to_BEDD
            
    allRRP = plant.getOrganRandomParameter(2)


    Main_beta = 3.28 * 10e-5
    Sub_beta = 2.16 * 10e-5
    SubSub_beta = 2.00 * 10e-5

    p2 = plant.getOrganRandomParameter(2, 2)
    allRRP[2].r = allRRP[2].lmax * Main_beta

    allRRP[2].tropismT = 7  # mix of planar and gravitropism
    allRRP[2].tropismN = 1  
    allRRP[2].tropismW1 = 0.85 #0.4 # gravitropism
    allRRP[2].tropismW2 = 0.15 #0.6 # plagiotropism

    allRRP[3].r = allRRP[3].lmax * Sub_beta
    
    for newpId in range(3,6):
        rrpA = allRRP[newpId].copy(plant)
        
        rrpA.subType += 1
        
        # will need to change that, to get the lmax ad so on from Lukas
        if rrpA.subType == 4:
            rrpA.k_survive = 2.17
            rrpA.lambda_survive = 4.55
            rrpA.r = rrpA.lmax * SubSub_beta
            
        rrpA.successorST =  np.array(rrpA.successorST) + 1
        allRRP.append(rrpA)
        

    allRRP[6].successorST = [[11]]
    allRRP[6].successorP = [[allRRP[6].successorP[0][1]]] # [[pp1],[pp2]]

    #pplats = np.array([[0.5]])#,[0.35]])
    for newpId in range(6,11):
        rrpA = allRRP[newpId].copy(plant)
        
        rrpA.subType += 1
        snext = rrpA.subType + 1
        rrpA.successorST =  [[snext]] # [[snext],[snext]]
        #pp = rrpA.successorP[0][0] * 0.5
        #rrpA.successorP =  pplats
        #rrpB.successorP =  np.array(rrpB.successorP) /10.
        rrpA.is_fine_root = True
        rrpA.lmax = 5. #ps.Lmax_suberized
        allRRP.append(rrpA)
        
        #pplats *= 0.75

    allRRP[-1].successorST = []

    for newp in allRRP:
        #newp.tropismS = 0.8
        plant.setOrganRandomParameter(newp)
       

        
        
    '''
    start simulation
    '''
    # TODO: select randomly one of the rsml files
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
      

    all_lengths = []
    all_ages = []
    all_subtypes = []
    all_alive = []
    dt = yr_to_BEDD  # ~1 yr
    N = 49

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
        
            
        # orrs = plant.getOrgans(-1)[0]
        # print(orrs.getId(),orrs.organType(),orrs.param().subType,orrs.getAge(),orrs.getLength(),
                    # orrs.getRootRandomParameter().a_gr, orrs.getRadius(0),
                    # # orrs.getNodeCT(0))
        # if i == 2:
            # raise Exception
        
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

