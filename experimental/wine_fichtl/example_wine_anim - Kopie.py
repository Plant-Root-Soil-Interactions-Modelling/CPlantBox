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

class My_SubType_Tropism(pb.Tropism):
	def __init__(self, plant, n, sigma):
		super(My_SubType_Tropism, self).__init__(plant)
		self.plagio =  pb.Plagiotropism(plant, n, sigma) 
		self.gravi = pb.Gravitropism(plant, n, sigma)
		self.setTropismParameter(n, sigma)
	def tropismObjective(self, pos, old, a, b, dx, root):
		subType = root.param().subType 
		if subType == 2:
			if np.random.rand() < 0.60:
				return self.gravi.tropismObjective(pos, old, a, b, dx, root)
			else:
				return self.plagio.tropismObjective(pos, old, a, b, dx, root)
				return self.gravi.tropismObjective(pos, old, a, b, dx, root)

plant = PlantPython(1)

# Open plant and root parameter from a file
path = "../../modelparameter/structural/rootsystem/"
name = "wine_Fichtl"

plant.readParameters(path + name + ".xml")

Main_beta = 3.28 * 10e-5
Sub_beta = 2.16 * 10e-5
SubSub_beta = 2.00 * 10e-5

#p to 3 a: 𝑋=0.35, 𝑌=0.55, 𝑍=0.10
#> 3 a: 𝑋=0.45, 𝑌=0.35, 𝑍=0.20; right?

ps = plant.getOrganRandomParameter(pb.seed)[0]
ps.Lmax_unsuberized = 5.
ps.Lmax_suberized = 10.
ps.delayDefinition = 4

p1 = plant.getOrganRandomParameter(2, 1)
p1.r = p1.lmax * Main_beta
p3 = plant.getOrganRandomParameter(2, 3)
p3.r = p3.lmax * SubSub_beta

p2 = plant.getOrganRandomParameter(2, 2)
#p2.successor = [[3]]
#p2.successorP = [[1]]
#print(p2.successorP, p2.successor, p2.a_gr, p2.r, p2.a)
#aise Exception
p2.r = p2.lmax * Sub_beta
#p2.ldelay = (300. * 60. )/2.
#p2.ldelays =(300. * 60. * 2.)/2.

# p2.theta = 0
# p2.thetas = 0
p2.tropismT = 1  #  1 gravi, 2 exo
p2.tropismS = 0.2
p2.tropismN = 0.5  # 0.05
p2.tropismW1 = 0.4
p2.tropismW2 = 0.6
p2.a_gr = p2.a/(300*60)
p2.is_fine_root = 0

p3 = plant.getOrganRandomParameter(2, 3)
p3.is_fine_root = 1

for ii, pp in enumerate(plant.getOrganRandomParameter(2)):
    if ii < 2:
        pp.ldelay  = 0*1e6
        pp.ldelays = 1e6*50
    else:
        pp.ldelay  = 0*1e6
        pp.ldelays = 1e6*50#/2.
allRRP = plant.getOrganRandomParameter(2)#[:2]
allRRP[2].successorST =  np.array(allRRP[2].successorST) + 2
allRRP[3].successorST =  np.array(allRRP[3].successorST) + 2
for newpId in range(2,10,2):
    rrpA = allRRP[newpId].copy(plant)
    rrpB = allRRP[newpId + 1].copy(plant)
    #print('init', rrpA.subType, rrpB.subType)
    rrpA.subType += 2
    rrpB.subType += 2
    rrpA.successorST =  np.array(rrpA.successorST) + 2
    rrpB.successorST =  np.array(rrpB.successorST) + 2
    rrpA.successorP =  np.array(rrpA.successorP) /2.
    rrpB.successorP =  np.array(rrpB.successorP) /10.
    allRRP.append(rrpA)
    allRRP.append(rrpB)
    #print('after', rrpA.subType, rrpB.subType)

# print(allRRP[-1].subType)
# print(allRRP[-2].subType)
for newp in allRRP:
    plant.setOrganRandomParameter(newp)
# for ii, pp in enumerate(plant.getOrganRandomParameter(2)):
    # print(pp.subType)
    
allRRP[-1].successorST = []
allRRP[-2].successorST = []

for newp in allRRP:
    plant.setOrganRandomParameter(newp)
 
for ii, pp in enumerate(plant.getOrganRandomParameter(2)):
    print(pp)
    
 
plant.initialize_static(path + "B-23_Fichtl.rsml", [0, 1])  # 0 is shoot, 1 are static roots

# the static laterals 2 and 3 are replaced with the growing lateral 2
plant.set_identical_laterals([0, 1], [2, 3], 2)
plant.initialize_static_laterals()

#p2.f_tf = My_SubType_Tropism

# plant.simulate(125., True)
dt = 1e6 # ~1 yr
N = 5
min_ = np.array([-20, -20, -50])
max_ = np.array([20, 20, 30.])

# test = plant.getOrgans(pb.leaf)
# print("test")

# anim = vp.AnimateRoots(plant)
# anim.min = min_
# anim.max = max_
# anim.res = [1, 1, 1]
# anim.file = "results/example_plant"
# anim.avi_name = "results/example_"
# anim.plant = True
# anim.start()
rlt_winter = plant.getParameter('rlt_winter')
print('rlt_winter',rlt_winter)
indexI = 1
for i in range(0, N):

    plant.simulate(dt, False)
    plant.survivalTest()
    
ana = pb.SegmentAnalyser(plant) 
vp.plot_roots(plant, "subType")
sts = np.array(plant.getParameter('subType', all = True))
creationTime = np.array(plant.getParameter('creationTime', all = True))
print(np.where([st in [3,5,7,9,11] for st in sts])[0])
plt.figure()
plt.hist([creationTime[np.where([st in [3,5,7,9,11] for st in sts])[0]],
            creationTime[np.where([st in [2,4,6,8,10] for st in sts])[0]]], stacked=True, density = False, bins=30, edgecolor='black',
            label=['fine roots', 'longlived roots'] )
#plt.hist(creationTime[sts == 3], stacked=True, density = True, bins=30, edgecolor='black')
plt.legend()
plt.show()
# print(dt*N, len(plant.getOrgans(2,False)), len(plant.getOrgans(2,True)))
# raise Exception

# rlt_winter = plant.getParameter('rlt_winter')
# orgId = np.array(plant.getParameter('id'))
# #print('rlt_winter',rlt_winter)
# #print('orgId',orgId)

# # anim.root_name = "organType"
# # anim.update()

# #plant.survivalTest()
# alives = plant.getParameter('alive')
# age = plant.getParameter('age')
# is_fine_root = plant.getParameter('is_fine_root')

# #print('alives',alives)
# #print('rlt_winter',rlt_winter)
# #print('age',age)
# #print('is_fine_root',is_fine_root)
# alives = plant.getParameter('alive')

# ana = pb.SegmentAnalyser(plant) 
# ana.addData('alive',alives)
# vp.plot_roots(plant, "alive")
# vp.plot_roots(plant, "subType")
