from fipy import *
from fipy.meshes.nonUniformGrid1D import NonUniformGrid1D as Grid1D
from fipy.meshes.nonUniformGrid2D import NonUniformGrid2D as Grid2D
from fipy.tools import numerix as np
import sys; sys.path.append("../../.."); sys.path.append("../../../src/python_modules")
from xylem_flux import XylemFluxPython  # Python hybrid solver
from phloem_flux import PhloemFluxPython, GetTime
from phloem_fluxbis import PhloemFluxPythonbis

from CellVariablemod import CellVariablemod
import plantbox as pb
import vtk_plot as vp
import math
import os
from io import StringIO
from datetime import datetime, timedelta

home_dir = os.getcwd()
dir_name = "/results"
dir_name2 = home_dir + dir_name
test = os.listdir(dir_name2)
for item in test:
    if item.endswith("9a.txt"):
        os.remove(os.path.join(dir_name2, item))
    if item.endswith("9a.vtk"):
        os.remove(os.path.join(dir_name2, item))
    if item.endswith("9a.vtp"):
        os.remove(os.path.join(dir_name2, item))
        
        
np.set_printoptions(threshold=sys.maxsize)

class NullIO(StringIO):
    def write(self, txt):
       pass

######################
#
# plant
#
####################### 
pl = pb.MappedPlant() #pb.MappedRootSystem() #pb.MappedPlant()

path = "../../../modelparameter/plant/" #"../../../modelparameter/rootsystem/" 
name = "smallPlant_mgiraud"#"manyleaves"#"oneroot_mgiraud" #"manyleaves"
pl.readParameters(path + name + ".xml")
start =1

""" Parameters xylem"""
kz = 4.32e-1  # axial conductivity [cm^3/day] 
kr = 1.728e-4  # radial conductivity of roots [1/day]
kr_stem = 1.e-20  # radial conductivity of stem  [1/day], set to almost 0
gmax = 0.004 #  cm3/day radial conductivity between xylem and guard cell
p_s = -200  # static water potential (saturation) 33kPa in cm
#p_g = -2000 # water potential of the guard cell
RH = 0.5 # relative humidity
TairC = 20
p_a =  -1000  #default outer water potential 
simtime = 15 # [day] for task b
k_soil = []
Q = 900e-6 # mol quanta m-2 s-1 light, example from leuning1995
cs = 350e-6 #co2 paartial pressure at leaf surface (mol mol-1)
TairK = TairC + 273.15


es = 0.61078 * math.exp(17.27 * TairC / (TairC + 237.3)) 
ea = es * RH 
VPD = es - ea 


""" Parameters phloem """
""" soil """
min_ = np.array([-50, -50, -150])
max_ = np.array([90, 40, 10])
res_ = np.array([5, 5, 5])
pl.setRectangularGrid(pb.Vector3d(min_), pb.Vector3d(max_), pb.Vector3d(res_), True)  # cut and map segments
pl.initialize()

pl.simulate(start, False)
phl = PhloemFluxPython(pl)
segs = pl.getPolylines() #segments regrouped per organ
phl.mp2mesh(segs) #creates grid
organTypes = phl.get_organ_types()
print(phl.mesh.length, "\n", phl.orgID, "\n", phl.orgLength, "\n", [organTypes[phl.new2oldNodeID[xi] - 1] for xi in phl.mesh.cellFaceIDs[1]])

print("2nd sim")
pl.simulate(start, False)
print("\n2nd sim end")
phl = PhloemFluxPython(pl)
phlbis = PhloemFluxPythonbis(pl)
segs = pl.getPolylines() #segments regrouped per organ
phl.mp2mesh(segs) #creates grid
organTypes = phl.get_organ_types()
print(phl.mesh.length, "\n", phl.orgID, "\n", phl.orgLength, "\n", [organTypes[phl.new2oldNodeID[xi] - 1] for xi in phl.mesh.cellFaceIDs[1]])


print("3rd sim")
pl.simulate(start, False)
print("\n3rd sim end")
phl = PhloemFluxPython(pl)
phlbis = PhloemFluxPythonbis(pl)
segs = pl.getPolylines() #segments regrouped per organ
phl.mp2mesh(segs) #creates grid
organTypes = phl.get_organ_types()
print(phl.mesh.length, "\n", phl.orgID, "\n", phl.orgLength, "\n", [organTypes[phl.new2oldNodeID[xi] - 1] for xi in phl.mesh.cellFaceIDs[1]])
print(phl.mesh.cellFaceIDs)

cellsIDOld = [phl.new2oldNodeID[xi] - 1 for xi in phl.mesh.cellFaceIDs[1] ]


######################

#tip indexes 
tiproots, tipstem, tipleaf = phl.get_organ_nodes_tips()
tiproots_newID = np.array([phl.old2newNodeID[xi] for xi in tiproots]).flatten()
tiproots=  [np.any([tiproots_newID == yi]) for yi in phl.mesh.faceVertexIDs[0]]
tiprootsfaces= sum([phl.mesh.faceVertexIDs[0] == xi for xi in tiproots_newID]) 
tiprootsfacesID= np.where(tiprootsfaces)[0]
tiprootscells= sum([phl.mesh.cellFaceIDs[1] == xi for xi in tiprootsfaces])
tiprootscellsID= np.where(tiprootscells)[0]

####
#
# Equations
#
####


cumulOutOld = CellVariable(mesh = phl.mesh, value=0.)
cumulGrOld = CellVariable(mesh = phl.mesh, value=0.)
cumulRmOld = CellVariable(mesh = phl.mesh, value=0.)
cumulAnOld = CellVariable(mesh = phl.mesh, value=0.)
phl.resetValues()    
dt =  1#0.9 * min(phl.mesh.length) ** 2 / 2#(2 * max(phl.Source)) #set to seconds
phiOld = CellVariablemod(mesh = phl.mesh, value=0., hasOld = True)#180 * dt*0.1

cellVolumeOld = phl.mesh.cellVolumes

logfilermM = open('results/rmMax_9a.txt', "w")
logfilerm = open('results/rm_9a.txt', "w")
logfilephi = open('results/phi_9a.txt', "w")
logfilerg = open('results/rg_9a.txt', "w")
logfilergS = open('results/rgSink_9a.txt', "w")
logfileOut = open('results/outFow_9a.txt', "w")


phl.phi.updateOld()
growthSteps = []
issue = []
issueRes = []
issueLoop = []
phiConcentrationError = []


phl.oldCellsID = [phl.new2oldNodeID[xi] - 1 for xi in phl.mesh.cellFaceIDs[1] ]


timeSpanLeuning = 1 #every 30mn #time in seconds between to evaluation of the gs operture/xylem flow
timeSinceLastLeuning = 0#np.inf

timeSinceLastGrowth = 0#np.inf
organGr = np.array([0.])

#for _ in range(steps):
step = 0
beginning = datetime.now()
simDuration = 0


while simDuration < 20:#60*60*24: #
    print('step n°', step,', dt: ',dt,'s, ','tot sim time: ', GetTime(simDuration))

    #sys.stdout = NullIO()  
  
    #step n° 4776 , dt:  10 s,  tot sim time:  0d 13hr 15mn 51s  : error
    if (timeSinceLastLeuning > timeSpanLeuning or step == 0):
          ####
        #
        # plant growth
        #
        ####
        print('growth' , phl.orgID)
        if(len(organGr) > 1):
            pl.setCWGr(organGr ) 
        pl.simulate(timeSinceLastGrowth/(60*60*24) , False) #time span in days /(60*60*24)
        
        phl = PhloemFluxPython(pl)
        segs = pl.getPolylines() #segments regrouped per organ
        timeSinceLastGrowth  = 0
        print('growthEnd')
        ####
        #
        # xylem
        #
        ####
        
        nodes = phl.get_nodes()
        tiproots, tipstem, tipleaf = phl.get_organ_nodes_tips() #end node of end segment of each organ
        node_tips = np.concatenate((tiproots, tipstem, tipleaf))
        tiproots, tipstem, tipleaf = phl.get_organ_segments_tips() #end segment of each organ
        seg_tips = np.concatenate((tiproots, tipstem, tipleaf))


        phl.setKr([[kr],[kr_stem],[gmax]]) #att: check units
        phl.setKx([[kz]])
        phl.airPressure = p_a

        phl.seg_ind = seg_tips # segment indices for Neumann b.c.
        phl.node_ind = node_tips
        rx = phl.solve_leuning( sim_time = dt,sxx=[p_s], cells = True, Qlight = Q,VPD = VPD,
            Tl = TairK,p_linit = p_s,ci_init = cs,cs=cs, soil_k = [], log = False)
        fluxes = phl.radial_fluxes(timeSinceLastLeuning/(60*60*24), rx, [p_s], k_soil, True)  # cm3/day
        
        timeSinceLastLeuning  = 0
        organGr = [0]
        
        phl.Px = rx
    
        ####
        #
        # phloem
        #
        ####
        print('mesh before ')
        phl.mp2mesh(segs) #creates grid
        print('mesh end ')
        phl.cellsID = [phl.new2oldNodeID[xi] - 1 for xi in phl.mesh.cellFaceIDs[1] ]
        print(phl.mesh.cellFaceIDs, phl.mesh.faceVertexIDs,  "\n", phl.orgLength, phl.mesh.length)
        if(len(phl.phi.value) > len(phiOld)): #place according to index ID
            print('bigger ', len(phl.phi.value) , len(phiOld.value))
            CphiOld = phiOld * cellVolumeOld
            print('cellsID ',cellsIDOld,phl.cellsID)
            orderedCPhi = np.concatenate((np.array([x for _,x in sorted(zip(cellsIDOld, CphiOld))]), np.full(len(phl.phi.value)-len(phiOld.value),0.)))        
            orderedCPhiNew = np.take(orderedCPhi, phl.cellsID)/phl.mesh.cellVolumes
            growthSteps = np.append(growthSteps, [step])
            phl.phi = CellVariablemod(mesh = phl.mesh, value= orderedCPhiNew , hasOld =True)
            
            phl.phi.updateOld()
            
            orderedcumulOut = np.concatenate((np.array([x for _,x in sorted(zip(cellsIDOld, cumulOutOld))]), np.full(len(phl.phi.value)-len(phiOld),0.)))        
            orderedcumulOutNew = np.take(orderedcumulOut, phl.cellsID)        
            cumulOut = CellVariable(mesh = phl.mesh, value= orderedcumulOutNew)
            
            orderedcumulGr = np.concatenate((np.array([x for _,x in sorted(zip(cellsIDOld, cumulGrOld))]), np.full(len(phl.phi.value)-len(phiOld),0.)))        
            orderedcumulGrNew = np.take(orderedcumulGr, phl.cellsID)                
            cumulGr = CellVariable(mesh = phl.mesh, value= orderedcumulGrNew)
            
            orderedcumulRm = np.concatenate((np.array([x for _,x in sorted(zip(cellsIDOld, cumulRmOld))]), np.full(len(phl.phi.value)-len(phiOld),0.)))        
            orderedcumulRmNew = np.take(orderedcumulRm, phl.cellsID)                
            cumulRm = CellVariable(mesh = phl.mesh, value= orderedcumulRmNew)
            
            orderedcumulAn = np.concatenate((np.array([x for _,x in sorted(zip(cellsIDOld, cumulAnOld))]), np.full(len(phl.phi.value)-len(phiOld),0.)))        
            orderedcumulAnNew = np.take(orderedcumulAn, phl.cellsID)                
            cumulAn = CellVariable(mesh = phl.mesh, value= orderedcumulAnNew)
        
        else: 
            phl.phi = CellVariablemod(mesh = phl.mesh, value= phiOld.value * cellVolumeOld/phl.mesh.cellVolumes, hasOld =True)
            cumulOut = cumulOutOld
            cumulGr =cumulGrOld
            cumulRm = cumulRmOld
            cumulAn = cumulAnOld
    print('reset')        
    phl.resetValues()    
    
    print('eq')
    #um of ink term hold tay <= 0.8 of phl.phi. otherwise, issue with convergence
    #in work of uys, diff term on other side...
    eq = (TransientTerm(var = phl.phi) ==  (phl.phi.faceValue * phl.phi.faceGrad*phl.intCoeff ).divergence + phl.Source + ImplicitSourceTerm(var = phl.phi,coeff=-phl.totalSinkFactor/dt )-phl.outFlow)
        #- phl.totalSinkFactor )
        #
        #+ ImplicitSourceTerm(var = phl.phi,coeff=- phl.StructC* phl.k1))# - phl.StructC * phl.k2)
    eq2 = (TransientTerm(var = phl.Rm) ==  -ImplicitSourceTerm(var = phl.phi,coeff=-phl.totalSinkFactor/dt ))
    #eq2= (TransientTerm(var = phl.Rm) ==  TransientTerm(var = phl.phi) - phl.Source)
    ### -phl.Rm)# - phl.outFlow- phl.GrSink ) + ImplicitSourceTerm(var = phl.phi,coeff=-phl.totalSinkFactor )
    eqt = eq & eq2
    res = 1e+10
    resOld = 2e+10
    loop = 0
    
    simDuration += dt
    print('eqend')
    
    #phibeforesweep = phl.phi.value
    while ( (res > 1e-6*max(phl.phi* phl.mesh.cellVolumes)  or any(phl.phi <0)) and loop < 100 ) : #and resOld != res
        resOld = res
        res = eqt.sweep(dt= dt) #quid solver?
        #FiPy’s order of precedence when choosing the solver suite for generic solvers is Pysparse followed by Trilinos, PyAMG and SciPy.
        #and PyAMGX?
        #solver option with scipy:
        #DefaultSolver = LinearLUSolver
        #DummySolver = LinearGMRESSolver
        #DefaultAsymmetricSolver = LinearLUSolver
        #GeneralSolver = LinearLUSolver
        #__all__.extend(linearCGSSolver.__all__) => LinearCGSSolver => 3s
        #__all__.extend(linearGMRESSolver.__all__) => LinearGMRESSolver =>  3
        #__all__.extend(linearBicgstabSolver.__all__) => not defined
        #__all__.extend(linearLUSolver.__all__) => 3
        #__all__.extend(linearPCGSolver.__all__) => 3s
        #help for solver: https://www.mail-archive.com/fipy@nist.gov/msg01089.html
      
        loop += 1
        print('res ',res)#, phl.RmMax/phl.phi <= 1, phl.Rm)#, phl.phi, phl.Rm)#, phl.Rm, phl.GrSink,phl.outFlow , phl.Source)
    
    #print('phi ',  phl.phi, phl.Rm)
    if res > 1e-6*max(phl.phi* phl.mesh.cellVolumes) :
        print('no convergence! res: ', res)
        issue = np.append(issue, [step])
        issueRes = np.append(issueRes, [res])
        issueLoop = np.append(issueLoop, [loop])
    phl.Rm.updateOld()
    phl.phi.updateOld()
    #phl.Rm = phibeforesweep +  
    #print(phl.phi , phl.totalSinkFactor, phl.totalSinkFactor/(phl.RmMax + phl.Gr) ,  phl.totalSinkFactor*phl.phi, phl.RmMax , phl.Gr)
    #eq2.solve(dt = dt)
    #print('Rm ', phl.Rm)
    #cumulRm.setValue(cumulRm.value +phl.totalSinkFactor * dt* phl.mesh.cellVolumes)
    #cumulRm.setValue(cumulRm.value + phl.phi*(phl.totalSinkFactor/dt) * dt* phl.mesh.cellVolumes)#(cumulRm.value + phl.phi*phl.Rm.value * dt* phl.mesh.cellVolumes)
    cumulRm.setValue(cumulRm.value +phl.Rm.value * phl.mesh.cellVolumes)
    cumulAn.setValue(cumulAn.value + phl.Source.value * dt* phl.mesh.cellVolumes)
    #cumulGr.setValue(cumulGr.value + phl.GrSink.value * dt* phl.mesh.cellVolumes)
    cumulOut.setValue(cumulOut.value + phl.outFlow.value * dt* phl.mesh.cellVolumes)
    timeSinceLastLeuning += dt
    timeSinceLastGrowth += dt
    
    if step *dt % (60) == 0: #choose how often get output       
        ####
        #
        # copy of phi
        #
        ####
        
        reorderedPhi = np.array([x for _,x in sorted(zip(phl.cellsID,phl.phi.value))])
        reorderedOutFlow = np.array([x for _,x in sorted(zip(phl.cellsID,cumulOut))])
        reorderedRm = np.array([x for _,x in sorted(zip(phl.cellsID,cumulRm))])
        reorderedGrSink = np.array([x for _,x in sorted(zip(phl.cellsID,cumulGr))])
        reorderedVolume = np.array([x for _,x in sorted(zip(phl.cellsID,phl.mesh.cellVolumes))])
        ana = pb.SegmentAnalyser(phl.rs)
        ana.addData("phi",  np.around(reorderedPhi*1e10, 5)) 
        #need to do around or too many digit for paraiew. and add 1e10 otherwise just yields 0
        ana.addData("outFlow",  np.around(reorderedOutFlow/reorderedVolume*1e10, 5))
        ana.addData("Rm",  np.around(reorderedRm/reorderedVolume*1e10, 5))
        ana.addData("GrSink",  np.around(reorderedGrSink/reorderedVolume*1e10, 5))
        ana.addData("rx", phl.Px)
        ana.addData("fluxes", fluxes) 
        ana.write("results/%s_example9a.vtp" %(step), ["radius", "surface", "phi", "outFlow", "Rm", "GrSink", "rx","fluxes"])
        #print('phi4ana ',reorderedPhi, np.around(reorderedPhi*1e10, 5))
    if step  % 1 == 0: #choose how often get output    
        ####
        #
        # check balance
        #
        ####
        
        phimin = phl.phi * (phl.phi < 0 )
        #print('sink priority check ',sum(phl.GrSink * (~phl.CSat)),sum(phl.outFlow*~phl.satisfaction)) #tip of stem negative at the beginning
        print('sucrose concentration < 0 at ',phl.phi[np.where(phl.phi < 0 )[0]],np.where(phl.phi < 0 )[0], ' sucrose vol. fraction > 0.65 at ',
            len(np.where(phl.VolFractSucrose > 0.65 )[0]),phl.VolFractSucrose[np.where(phl.VolFractSucrose > 0.65 )[0]] )

        print('cumulated loading: ', sum(cumulAn), 'cumulated growth sink: ',sum(cumulGr), 'cumulated maintenance sink: ',sum(cumulRm), 'cumulated outflow: ',sum(cumulOut),'total C content in plant: ',sum(phl.phi* phl.mesh.cellVolumes))
        print('C_plant + C_Sink - C_source :',-sum(cumulAn) + sum(phl.phi* phl.mesh.cellVolumes)+ sum(cumulRm)+ sum(cumulOut) )# + sum(cumulGr)+ sum(cumulOut))
        diff = -sum(cumulAn) + sum(phl.phi* phl.mesh.cellVolumes)+ sum(cumulRm) + sum(cumulOut) + sum(cumulGr)
        #if diff/ sum(phl.phi* phl.mesh.cellVolumes) > 0.001:
         #   raise Exception("diff/ sum(phl.phi* phl.mesh.cellVolumes) > 0.001")
        #print('check if value updates ',sum(phl.CSat), sum(phl.RmMax), sum(phl.Rm), sum(phl.Gr), sum(phl.GrSink), sum(phl.outFlow)) #check if value updates
        #print( 'Growth sink in segments with C limitation for maintenance (should be 0): ',) #check if value updates
        '''
        logfilephi.write('\n'+repr(phl.phi.value)[7:-2])
        logfilerm.write('\n'+repr(phl.Rm.value/phl.mesh.cellVolumes)[7:-2])
        logfilermM.write('\n'+repr(phl.RmMax.value/phl.mesh.cellVolumes)[7:-2])
        logfilergS.write('\n'+repr(phl.GrSink.value/phl.mesh.cellVolumes)[7:-2])
        logfilerg.write('\n'+repr(phl.Gr/phl.mesh.cellVolumes)[7:-2])
        logfileOut.write('\n'+repr(phl.outFlow.value/phl.mesh.cellVolumes)[7:-2])
        '''
        print('no convergence at step(s) ',issue, ' res: ', issueRes, ' loops ', issueLoop)
        print('phi value error at step(s) ',phiConcentrationError)
        #print('became bigger at step(s) ',growthSteps)
    phiOld = phl.phi.copy()
    cumulOutOld = cumulOut.copy()
    cumulGrOld = cumulGr.copy()
    cumulRmOld = cumulRm.copy()
    cumulAnOld = cumulAn.copy()
    cellVolumeOld = phl.mesh.cellVolumes
    cellsIDOld = phl.cellsID
    print('growth per segment :',organGr ,'\ndt :', dt, phl.orgID)
    organGr += phl.organGr * dt
    print('growth per segment :',organGr ,'\ndt :', dt, phl.orgID)
    step += 1
    
    dt =  10#min((0.9 * min(phl.mesh.length) ** 2 / (2 * max((phl.phi.faceValue*phl.intCoeff).value)),30))
    #print('calc dt ',dt,"\n", phl.phi, phl.RmMax, phl.Rm, phl.Gr,phl.GrSink, phl.CSat, phl.outFlow ) #[-41:]
    #print(phl.mesh.cellFaceIDs, phl.mesh.faceVertexIDs, phl.mesh.length, phl.mesh.cellVolumes)
    if(any(phl.phi < 0)):
        raise Exception("phl.phi < 0")
    if(any([math.isnan(ci) for ci in phl.phi])):
        raise Exception("phl.phi is nan")
    if(any(phl.Rm < 0)):
        raise Exception("phl.phl.totalR < 0")
    
    #if(sum(abs(phl.GrSink * (~phl.CSat)))+sum(abs(phl.outFlow*~phl.satisfaction))\
     #   + sum(abs(phl.phi[np.where(phl.phi < 0 )[0]]))+sum(abs(phl.VolFractSucrose[np.where(phl.VolFractSucrose> 0.65 )[0]])) != 0):
      #  raise Exception("error priority or sucVolFraction")
    
end = datetime.now()
duration = end - beginning
print('duration simulation ', duration.seconds, 's for a simulation of ', step ,' steps and ', GetTime(simDuration),  "\n", phl.orgLength, phl.mesh.length)

logfilermM.close()
logfilerm.close()
logfilephi.close()
logfilerg.close()
logfilergS.close()
logfileOut.close()
