from fipy import *
from fipy.meshes.nonUniformGrid1D import NonUniformGrid1D as Grid1D
from fipy.meshes.nonUniformGrid2D import NonUniformGrid2D as Grid2D
from fipy.tools import numerix as np
import sys; sys.path.append("../../.."); sys.path.append("../../../src/python_modules")
from xylem_flux import XylemFluxPython  # Python hybrid solver
from phloem_flux import PhloemFluxPython, GetTime

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
    if item.endswith("10f.txt"):
        os.remove(os.path.join(dir_name2, item))
    if item.endswith("10f.vtk"):
        os.remove(os.path.join(dir_name2, item))
    if item.endswith("10f.vtp"):
        os.remove(os.path.join(dir_name2, item))
        
        
np.set_printoptions(threshold=sys.maxsize)

class NullIO(StringIO):
    def write(self, txt):
       pass
       

####
#
# wheather data
#
####
def weatherData(t):
    t_ = round(t)
    tmp = ((t_+ 12) % 24) + (t/(60*60) - t_)
    trapeze = max(0, min(1, min(tmp-11, 24-tmp))) 
    TairC =  20*trapeze+5
    TairK = TairC + 273.15
    RH = 0.5 + 0.4*(1-trapeze)# relative humidity
    Q = 900e-6*trapeze # mol quanta m-2 s-1 light, example from leuning1995
    
    

    es = 0.61078 * math.exp(17.27 * TairC / (TairC + 237.3)) 
    ea = es * RH 
    VPD = es - ea 

    weatherVar = {'TairK' : TairK,
                    'Qlight': Q,#(-)
                    'VPD': VPD}
    print('weather ', t,trapeze,  weatherVar)
    return weatherVar
    

    
    
######################
#
# plant
#
####################### 
pl = pb.MappedPlant() #pb.MappedRootSystem() #pb.MappedPlant()

path = "../../../modelparameter/plant/" #"../../../modelparameter/rootsystem/" 
name =  "Triticum_aestivum_adapted_2021"#"smallPlant_mgiraud"#"morning_glory_7m"#"manyleaves"#"oneroot_mgiraud" #"manyleaves"
pl.readParameters(path + name + ".xml")
start =35

""" Parameters xylem"""
kz = 4.32e-1  # axial conductivity [cm^3/day] 
kr = 1.728e-4  # radial conductivity of roots [1/day]
kr_stem = 1.e-20  # radial conductivity of stem  [1/day], set to almost 0
gmax = 0.004 #  cm3/day radial conductivity between xylem and guard cell
p_s = -200  # static water potential (saturation) 33kPa in cm
#p_g = -2000 # water potential of the guard cell
p_a =  -1000  #default outer water potential 
simtime = 15 # [day] for task b
k_soil = []
cs = 350e-6 #co2 paartial pressure at leaf surface (mol mol-1)
p_linit = -10000 #guess for first p_l value 

VariableMu = False


""" Parameters phloem """
""" soil """
min_ = np.array([-50, -50, -150])
max_ = np.array([90, 40, 10])
res_ = np.array([5, 5, 5])
pl.setRectangularGrid(pb.Vector3d(min_), pb.Vector3d(max_), pb.Vector3d(res_), True)  # cut and map segments
pl.initialize(True)

pl.simulate(start, False)


phl = PhloemFluxPython(pl)
ana = pb.SegmentAnalyser(phl.rs)
ana.write("results/example10f.vtp" )


       
phl.mp2mesh(VariableMu = VariableMu) #creates grid
organTypes = phl.get_organ_types()
cellsIDOld = [phl.new2oldNodeID[xi] - 1 for xi in phl.mesh.cellFaceIDs[1] ]
initphi =0.


initphi = 0.#np.loadtxt('Cfin.txt', skiprows=0)[1:]/1000/12 #mmol C to mol suc
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

phl.resetValues(0)    

phiOld = CellVariablemod(mesh = phl.mesh, value=  initphi, hasOld = True)#initial value (cf. xiaoran zhou) 
phiOld.updateOld()
SucSinkOld = phl.SucSink.copy()*0

MesophyllOld =phl.Mesophyll.copy()
MesophyllOld.setValue(initphi)
#MesophyllOld.updateOld()
cellVolumeOld = phl.mesh.cellVolumes

phl.phi.updateOld()
growthSteps = []
issue = []
issueRes = []
issueLoop = []
phiConcentrationError = []


phl.oldCellsID = [phl.new2oldNodeID[xi] - 1 for xi in phl.mesh.cellFaceIDs[1] ]


timeSpanLeuning = 5  #time in seconds between to evaluation of the gs operture/xylem flow
timeSinceLastLeuning = np.inf

timeSpanPlot = np.inf
timeSinceLastPlot =0.

timeSpanGrowth = 24 
timeSinceLastGrowth = 0.

timeSpanCheck = 1/6
timeSpanCheck = 1/6
timeSinceLastCheck = np.inf
numcheck = 0

organGr = [0.]

#for _ in range(steps):
step = 0
beginning = datetime.now()
simDuration = 0
dtVal = dt=0.125
alldts = np.array([])
maxdt = 60
increase = 1.5



        
while simDuration < 36: #
    end = datetime.now()
    duration = end - beginning
    print('\n\nstep n°', step,', dt: ',dt,'hrs, ','tot sim time: ', GetTime(simDuration), ', computation time: ',duration)
    #, 0.9 * min(phl.mesh.length)**2 / (2 * max(phl.phi.faceValue *phl.intCoeff)))

    
    if (timeSinceLastLeuning > timeSpanLeuning or step == 0):
        weatherD = weatherData(simDuration)
        #old_stdout = sys.stdout # backup current stdout
        #sys.stdout = open(os.devnull, "w")
          ####
        #
        # plant growth
        #
        ####
        #print('growth' , phl.orgID)
        
        #p_linit = mean(rx[])
        if (timeSinceLastGrowth >timeSpanGrowth):
            if(len(organGr) > 1):
                pl.setCWGr(organGr ) 
                organGr =[0.]
            pl.simulate(timeSinceLastGrowth/24 , False) #time span in days /(60*60*24)
            
            phl = PhloemFluxPython(pl)
            timeSinceLastGrowth  = 0
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
        rx = phl.solve_leuning( sim_time = dt/24,sxx=[p_s], cells = True, Qlight = weatherD['Qlight'] ,VPD = weatherD['VPD'],
            Tl = weatherD['TairK'] ,p_linit = p_linit,ci_init = cs*0.3,cs=cs, soil_k = [], log = False)
        fluxes = phl.radial_fluxes(timeSinceLastLeuning/24, rx, [p_s], k_soil, True)  # cm3/day
        
        timeSinceLastLeuning  = 0
        
        phl.Px = rx
    
        ####
        #
        # phloem
        #
        ####
        phl.mp2mesh(VariableMu = VariableMu) #creates grid
        phl.cellsID = [phl.new2oldNodeID[xi] - 1 for xi in phl.mesh.cellFaceIDs[1] ]
        
        if(len(phl.phi.value) > len(phiOld)): #place according to index ID
            #replace with function
            CphiOld = phiOld * cellVolumeOld
            orderedCPhi = np.concatenate((np.array([x for _,x in sorted(zip(cellsIDOld, CphiOld))]), np.full(len(phl.phi.value)-len(phiOld.value),0.)))        
            orderedCPhiNew = np.take(orderedCPhi, phl.cellsID)/phl.mesh.cellVolumes
            growthSteps = np.append(growthSteps, [step])
            phl.phi = CellVariablemod(mesh = phl.mesh, value= orderedCPhiNew , hasOld =True)
            
            phl.phi.updateOld()
            
            CSucSinkOld = SucSinkOld * cellVolumeOld
            orderedCSucSink = np.concatenate((np.array([x for _,x in sorted(zip(cellsIDOld, CSucSinkOld))]), np.full(len(phl.phi.value)-len(phiOld.value),0.)))        
            orderedCSucSinkNew = np.take(orderedCSucSink, phl.cellsID)/phl.mesh.cellVolumes
            phl.SucSink = CellVariablemod(mesh = phl.mesh, value= orderedCSucSinkNew , hasOld =True)
            phl.SucSink.updateOld()
            
            
            MesophyllOld = MesophyllOld * cellVolumeOld
            orderedMesophyll = np.concatenate((np.array([x for _,x in sorted(zip(cellsIDOld, MesophyllOld))]), np.full(len(phl.phi.value)-len(phiOld.value),0.)))        
            orderedMesophyllNew = np.take(orderedMesophyll, phl.cellsID)/phl.mesh.cellVolumes
            phl.Mesophyll = CellVariablemod(mesh = phl.mesh, value= orderedMesophyllkNew , hasOld =True)
            phl.Mesophyll.updateOld()
            
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
            phl.SucSink = CellVariablemod(mesh = phl.mesh, value= SucSinkOld.value * cellVolumeOld/phl.mesh.cellVolumes, hasOld =True)
            phl.Mesophyll = CellVariablemod(mesh = phl.mesh, value= MesophyllOld.value * cellVolumeOld/phl.mesh.cellVolumes, hasOld =True)
            cumulOut = cumulOutOld
            cumulGr =cumulGrOld
            cumulRm = cumulRmOld
            cumulAn = cumulAnOld
        #sys.stdout = old_stdout # reset old stdout   
    phl.resetValues(simDuration)    
    #print('check if value updates1 ',sum(phl.totalSinkFactor), sum(phl.RmMax),  sum(phl.outFlowMax)) #check if value updates
    #sum of sink term <= 0.8 of phl.phi. otherwise, issue with convergence
    
    convergence = False
    #print("js ",phl.Mesophyll,phl.Agross,  phl.JS_Mesophyll*dt,phl.Mesophyll+phl.JS_Mesophyll*dt-phl.Agross*dt)
    if(step > 0 ): #convergence criteria according to phi value before loop (as might get odd values during loop)
        phi4res = phl.phi.copy() - initphi
    else: #for first loop, phi == 0, so set phi4res as CellVariable => can update at each loop (otherwise convergence never reached)
        phi4res = phl.phi - initphi
    #print(phl.phloemConductivity *  np.pi*phl.mesh.radiiVertices**2 /phl.mesh._cellDistances/phl.phloemViscosity, phl.Param['R']*1e6*phl.TairK*( phl.phi.faceGrad.divergence)*12)
    #print('JW, ', (phl.phi.faceGrad*phl.intCoeff ).divergence)
    #print('JSdiff, ', (phl.phi.faceValue * phl.phi.faceGrad*phl.intCoeff ).divergence -phl.phi * (phl.phi.faceGrad*phl.intCoeff ).divergence)#does make a significant difference
    while (not convergence):
        eqMesophyll = (TransientTerm(var = phl.Mesophyll) == phl.Agross - phl.JS_Mesophyll)
        #eq = (TransientTerm(var = phl.phi) ==  (phl.phloemConductivity *  np.pi*phl.mesh.radiiVertices**2 /phl.mesh._cellDistances/phl.phloemViscosity*phl.Param['R']*1e6*phl.TairK*phl.phi.faceValue*12* phl.phi.faceGrad).divergence + phl.JS_Mesophyll -phl.totalSinkFactor)#+ ImplicitSourceTerm(var = phl.phi,coeff=-phl.totalSinkFactor ))
        #area and length already included in duvergence calculation
        #eq = (TransientTerm(var = phl.phi) ==  phl.phi * (phl.phi.faceGrad*phl.intCoeff ).divergence + phl.JS_Mesophyll -phl.totalSinkFactor)#+ ImplicitSourceTerm(var = phl.phi,coeff=-phl.totalSinkFactor ))
        #=> makes error
        eq = (TransientTerm(var = phl.phi) ==   (phl.phi.faceValue * phl.phi.faceGrad*phl.intCoeff ).divergence + phl.JS_Mesophyll -phl.totalSinkFactor)#+ ImplicitSourceTerm(var = phl.phi,coeff=-phl.totalSinkFactor ))
        #eq2 = (TransientTerm(var = phl.SucSink) ==  phl.totalSinkFactor)#ImplicitSourceTerm(var = phl.phi,coeff=phl.totalSinkFactor ))
        
        eqt = eq  & eqMesophyll #& eq2
        res =diff= 1e+10
        resOld = 2e+10
        loop = 0
        #or diff /max(phl.phi* phl.mesh.cellVolumes)> 0.1 
        while ( ( res > 1e-5*max(phi4res* phl.mesh.cellVolumes)or np.isnan(res) or any(phl.phi <0)or any(phl.SucSink <0)) and loop <= 100 ) : #and resOld != res, or (diff3 /max(phl.phi* phl.mesh.cellVolumes)> 0.0001 )
            resOld = res
            res = eqt.sweep(dt= dt) #quid solver?
            #diff = abs(sum(phl.phi* phl.mesh.cellVolumes) -sum( phiOld* phl.mesh.cellVolumes)-sum(self.Agross*dt* phl.mesh.cellVolumes)\
             #   + sum(phl.phi * phl.totalSinkFactor*dt* phl.mesh.cellVolumes))
            #print('res: ',res, ' lim = ', 1e-3*max(phi4res* phl.mesh.cellVolumes))
          
            loop += 1
            
        if loop > 100:#:res > 1e-6*max(phl.phi* phl.mesh.cellVolumes) :
            print('no convergence! res: ', res,' lim = ', 1e-5*max(phi4res* phl.mesh.cellVolumes))
            issue = np.append(issue, [step])
            issueRes = np.append(issueRes, [res])
            issueLoop = np.append(issueLoop, [loop])
            dtVal = dtVal/2
            dt = dtVal
            maxdt = 1
            increase = 1.01
            print("change dt from ", dtVal*2, " to ", dtVal)
        else:
            print('convergence reached')
            convergence = True
            alldts = np.append(alldts, dt)
                
    #print(sum(phl.JS_Mesophyll* phl.mesh.cellVolumes*dt),sum(phl.Mesophyll* phl.mesh.cellVolumes+phl.JS_Mesophyll*dt* phl.mesh.cellVolumes-phl.Agross*dt* phl.mesh.cellVolumes))
    #print('check1 ',sum(phl.SucSink), sum(phl.outFlowMax))
    #phl.SucSink.updateOld()
    phl.phi.updateOld()
    phl.Mesophyll.updateOld()
    print('jsmeso ', phl.JS_Mesophyll[np.where(phl.JS_Mesophyll>0)[0]]* phl.mesh.cellVolumes[np.where(phl.JS_Mesophyll>0)[0]])
    #,sum(phl.Mesophyll* phl.mesh.cellVolumes+phl.JS_Mesophyll*dt* phl.mesh.cellVolumes-phl.Agross*dt* phl.mesh.cellVolumes))
    phl.getCLimitedSinks(phl.totalSinkFactor.value, dt)# phl.SucSink.value - SucSinkOld.value, dt) #compute C-limited Rm, GrSink and outflow
    cumulRm.setValue(cumulRm.value +phl.Rm* dt * phl.mesh.cellVolumes)
    cumulAn.setValue(cumulAn.value + phl.JS_Mesophyll* dt* phl.mesh.cellVolumes)#phl.Agross.value * dt* phl.mesh.cellVolumes)
    cumulGr.setValue(cumulGr.value + phl.GrSink* dt * phl.mesh.cellVolumes)
    cumulOut.setValue(cumulOut.value + phl.outFlow* dt* phl.mesh.cellVolumes)
    phl.SucSink.setValue(phl.SucSink.value +phl.totalSinkFactor*dt*phl.mesh.cellVolumes)
    phl.SucSink.updateOld()
    
    simDuration += dt#; print('simDuration ',simDuration)
    timeSinceLastLeuning += dt
    timeSinceLastGrowth += dt
    timeSinceLastPlot += dt
    timeSinceLastCheck +=dt
    if timeSinceLastPlot >= timeSpanPlot: #choose how often get output       
        phl.printParaview(cumulAn, cumulOut,cumulRm,cumulGr, fluxes ,simDuration ) #vtp for paraview
        timeSinceLastPlot = 0.
        
    ####
    #
    # check balance
    #
    ####
    #print('check if value updates2 ',sum(phl.totalSinkFactor), sum(phl.RmMax),  sum(phl.outFlowMax)) #check if value updates
    phl.doCheck(cumulAn, cumulOut,cumulRm,cumulGr, initphi, timeSinceLastCheck , timeSpanCheck, numcheck,issue, issueRes,  issueLoop,phiConcentrationError, dt)
    
    if timeSinceLastCheck >= timeSpanCheck: #choose how often get output 
        numcheck +=1
        end = datetime.now()
        duration = end - beginning
        phl.doLog( duration, organTypes, dt, issue, issueRes, issueLoop, alldts, simDuration) #print in log
        timeSinceLastCheck = 0.
        
    phiOld = phl.phi.copy()
    SucSinkOld = phl.SucSink.copy()
    MesophyllOld = phl.Mesophyll.copy()
    cumulOutOld = cumulOut.copy()
    cumulGrOld = cumulGr.copy()
    cumulRmOld = cumulRm.copy()
    cumulAnOld = cumulAn.copy()
    cellVolumeOld = phl.mesh.cellVolumes
    cellsIDOld = phl.cellsID
    organGr = phl.organGr + np.concatenate((organGr, np.full(len(phl.organGr) - len(organGr),0.)))
    
    step += 1
    #if step ==10 :
     #   dtVal = 0.5
    
    if loop < 10:
        dtVal = min(dtVal*increase, maxdt)
    dt = dtVal#min(0.9 * min(phl.mesh.length)**2 / (2 * max(phl.phi.faceValue *phl.intCoeff)),1)# dtVal# 60*60#min((0.9 * min(phl.mesh.length) ** 2 / (2 * max((phl.phi.faceValue*phl.intCoeff).value)),30))
    
    phl.valueError() #raise exception if impossible values
    
    
end = datetime.now()
duration = end - beginning
print('duration simulation ', duration.seconds, 's for a simulation of ', step ,' steps and ', GetTime(simDuration),  "\n", dt)

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
