
import sys ;sys.path.append("../../.."); sys.path.append("../../../src/python_modules")
#sys.path.append("/home/rbtlm640/CPBW/CPlantBox/src/python_modules")
import numpy as np;

import plantbox as pb
import vtk_plot as vp

import importlib
importlib.reload(vp)
import matplotlib.pyplot as plt
from vtk_tools import *
import math
print(sys.path)
from Leuning import Leuning
import plantbox as pb
import vtk_plot as vp
import pandas as pd
from matplotlib.ticker import MaxNLocator
import numpy as np
import matplotlib.pyplot as plt
##parameters for example:
adaptSeed = True
adaptLeaf = True
adaptStem = True
leafRadial = True #radial or not
anim = False
zoomLeafShape = True
export = False
getImage = True


""" Parameters """
kz = 4.32e-1  # axial conductivity [cm^3/day] 
kr = 1.728e-4  # radial conductivity of roots [1/day]
kr_stem = 1.e-20  # radial conductivity of stem  [1/day], set to almost 0
gmax =  0.004 #  cm3/day radial conductivity of leaves = stomatal conductivity [1/day]
p_a =  -1000  #static air water potential 
simtime = 14.0  # [day] for task b
k_soil = []
plotResults = True
saveResults = False

t_init = 70
t_end = 90

# root system 
plant = pb.MappedPlant() #pb.MappedRootSystem() #pb.MappedPlant()
path = "../../../modelparameter/plant/" #"../../../modelparameter/rootsystem/" 
name = "manyleaves" #"Anagallis_femina_Leitner_2010"  # Zea_mays_1_Leitner_2010
plant.readParameters(path + name + ".xml")

#load data
df = pd.read_csv('../../../modelparameter/Selhausen_weather_data.txt', delimiter = "\t")

""" soil """
min_ = np.array([-5, -5, -15])
max_ = np.array([9, 4, 0])
res_ = np.array([5, 5, 5])
plant.setRectangularGrid(pb.Vector3d(min_), pb.Vector3d(max_), pb.Vector3d(res_), True)  # cut and map segments


# plant.initialize()
# plant.simulate(simtime, False)

soil_index = lambda x,y,z : max(int(np.floor(-z)),-1) #abovegroud nodes get index -1
plant.setSoilGrid(soil_index)

        

if adaptSeed:
    srp = pb.SeedRandomParameter(plant)  # with default values
    srp.firstTil = 0  # [day] first emergence of a tiler
    srp.delayTil = 0  # [day] delay between the emergence of tilers
    srp.maxTil = 0 # [-] number of tillers 
    plant.setOrganRandomParameter(srp)
    

if adaptStem:
    for p in plant.getOrganRandomParameter(pb.stem):
        if (p.subType > 0): # can be changed according to the suptypes of the plant
            p.nodalGrowth = 1   #< whether to implement the internodal growth 
            p.delayLat = 1  #< delay between stem creation and start of nodal growth [day]
            p.delayNG = 10   #< delay between lateral creation and growth [day]
            #p.tropismAge = 10 #< only used if tropsimT = 6
            plant.setOrganRandomParameter(p)
            
if adaptLeaf:
    for p in plant.getOrganRandomParameter(pb.leaf):
    
                #p.lmax - p.la - p.lb = leafMid = center of radial circle
        if (p.subType >= 2): #leaf subtypes start at 2
            p.lb =  1 # length of leaf stem
            p.la,  p.lmax = 3.5, 8.5
            p.areaMax = 10  # cm2, area reached when length = lmax
            N = 100  # N is rather high for testing
            if leafRadial:
            
                #LongLeaf:
                p.lb =  1 # length of leaf stem
                p.ln = 0
                p.la,  p.lmax = 3.5, 8.5
                p.areaMax = 10  # cm2, area reached when length = lmax
                N = 100  # N is rather high for testing
                phi = np.array([-90, -45, 0., 45, 90]) / 180. * np.pi
                l = np.array([3, 2.2, 1.7, 2, 3.5]) #distance from leaf center
                
                p.createLeafRadialGeometry(phi, l, N)
#                 print('')
                
            else:
                p.lb =  2 # length of leaf stem
                p.la,  p.lmax = 3.5, 8.5
                p.areaMax = 10  # cm2, area reached when length = lmax
                N = 100  # N is rather high for testing  #Semble non-fonctionnel
                y = np.array([-3, -3 * 0.7, 0., 3.5 * 0.7, 3.5])
                l = np.array([0., 2.2 * 0.7, 1.7, 1.8 * 0.7, 0.])
                p.createLeafGeometry(y, l, N)     
                
                
            p.tropismT = 6 # 6: Anti-gravitropism to gravitropism
            #p.tropismN = 5
            #p.tropismS = 0.1
            p.tropismAge = 10 #< age at which tropism switch occures, only used if p.tropismT = 6
            plant.setOrganRandomParameter(p)


for p in plant.getOrganRandomParameter(pb.leaf):
    if (p.subType >= 2):
        #print(p) #permet de print tous les paramètres de la feuille   
        plant.initialize()

if anim:
    dt = 1
    N_ = 50
    min_ = np.array([0, -20, 0])/2
    max_ = np.array([20, 20, 30.])/2
    anim = vp.AnimateRoots(plant)
    anim.min = min_
    anim.max = max_
    anim.res = [1, 1, 1]
    #anim.file = "results/example_plant"
    #anim.avi_name = "results/example_"
    anim.plant = True
    anim.start()
    for i in range(0, N_):
        plant.simulate(dt, False)
        anim.root_name = "subType"
        anim.update()

if getImage:
    # Simulate
    if not anim:
        plant.simulate(simtime, True)
        soil_index = lambda x,y,z : max(int(np.floor(-z)),-1) #abovegroud nodes get index -1
        plant.setSoilGrid(soil_index)
        
    
    #Plotting of the plant and possible creation of ".obj" files for further uses.
    outputDirectory = "output"
#     vp.plot_plant(plant, "organType") #Default
    NormalsZValue= []
    vp.plot_plant(plant, "organType", render = False, printout = True, outputDirectory = outputDirectory, timestamp = False, sim_name = "", GraphicalAccuracy = True, NormalsZValue = NormalsZValue)
#     print("NormalsZValue",len(NormalsZValue), NormalsZValue)
    
    
    #Oeuf de Pâques:
#     vp.plot_plant(plant, "organType", render = False, printout = True, sim_name = "Oeuf", date = "Pâques", outputDirectory = outputDirectory)
    
    #Creates the lists needed for further use of the computed light onto the plant.
    stems = plant.getOrgans(pb.stem)
    StemSegIDList = []
    for stem in stems:
        for i in range(0, stem.getNumberOfNodes() - 1):
            IdSeg = (stem.getNodeId(i+1)-1)  #allows the retrieval of the absolute node IDs from the organ.
            StemSegIDList.append(IdSeg)
#     print("StemSegIDList", len(StemSegIDList), StemSegIDList)
    
    leafes = plant.getOrgans(pb.leaf)
    NumberOfLeafes = len(leafes)
    print("Number of leaves", NumberOfLeafes)
    
    #Copy of what's done in the vtk_plot.py code in the definition "plot_plant" with LeafSegIDList as vp.LeafSegIDListList
    #As shown below, it can be either retrieved or remade easily.
    LeafSegIDList = vp.LeafSegIDListList
#     LeafSegIDList = []
#     for leaf in leafes:
#         for i in range(0, leaf.getNumberOfNodes() - 1):
#             IdSeg = (leaf.getNodeId(i+1)-1)  #allows the retrieval of the absolute node IDs from the organ.
#             LeafSegIDList.append(IdSeg)
#     print("LeafSegIDList", len(LeafSegIDList), LeafSegIDList)
    
    
    #Retrieve the thickness of the leaf. The leaf thickness is egal to the radius of the segments tagged "leaf".
#     print(np.array(plant.getParameter("radius"))) #Returns the entirity of the radiuses used in the organism.
    leafes = plant.getOrgans(pb.leaf)
    LeafThicknesses = []
    for l in leafes:
        Thickness =l.getParameter("radius")
        LeafThicknesses.append(Thickness)
#         print(radius)
#     print(LeafThicknesses)
    
    if not(outputDirectory.endswith("/")):
        outputDirectory = outputDirectory + "/"
    with open(outputDirectory + "/" + 'Thickness.txt', 'w') as f: #Stores the different leaf thicknesses.
        for i in range(len(LeafThicknesses)):
            f.write(str(LeafThicknesses[i]))
            if i != len(LeafThicknesses)-1: #writes to the next line and doesn't add a blank line at the bottom.
                f.write('\n') 
    f.close()
    

    #Creates a render window of the plant based on the amount of light received on each part via the light algorithm.   
    Analyser = pb.SegmentAnalyser(plant)
    SegsLength = len(np_convert(Analyser.segments))
    print("Number of segments: ", SegsLength)
    Test = []
#     Test = np.random.random((SegsLength,1))
    Test = np.zeros((SegsLength))
    StemNormalsZValue = np.zeros((SegsLength))
#     print("StemNormalsZValue", len(StemNormalsZValue), StemNormalsZValue)
    
    for i in range(Test.shape[0]):
        Test[i] = i
#         Test[i] = Test[i]/5000
    Test_Stem = Test
#     print("Test_Stem", len(Test_Stem),Test_Stem)
    
    QuadCounter = 0
    QuadIdList = []  #This list will contain the indices between which the read value from the input file (Irradiance.txt) must be written at (their are more than one quadrilateral per segment).
    #Note: the next line doesn't work if the code didn't go through a "plot_plant" first to create the "vp.QuadCounterListList"
    for c in vp.QuadCounterListList:
        QuadIdList.append(QuadCounter)
        QuadCounter = QuadCounter + c
    QuadIdList.append(QuadCounter)
#     print("QuadIdList",len(QuadIdList), QuadIdList)
    
    Test_Leaf = np.zeros((QuadCounter)) #QuadCounter is equal to the total number of quadrilaterals in the leaves.
    for i in range(len(QuadIdList)-1):
        for l in range(QuadIdList[i], QuadIdList[i+1]):
            Test_Leaf[l] = i #Simple increment for test purposes.
#     print("Test_Leaf",len(Test_Leaf),Test_Leaf)
    
    #Reading irradiance data from file
    Testing = True
    if Testing:
        Irradiance_Leaf = np.zeros((QuadCounter)) #creates empty numpy arrays of the right sizes to export in "ExtraParam" of plot_plant
        Irradiance_Stem = np.zeros((SegsLength))
        
        #the following code allows the user to retrieve data written after the light computation
#         import glob
#         if  not(outputDirectory.endswith("/")) and (outputDirectory != ""): #if the directory name is empty, it doesn't add a "/"
#             outputDirectory = outputDirectory + "/"
#         pathname = glob.glob("./" + outputDirectory + 'Irradiance.txt')
#     #     print("pathname:", pathname[0])
#         with open(pathname[0]) as f:
#             lines = f.read()
#             ValStemList = [0]*SegsLength
#             ValLeafList = [0]*(len(QuadIdList)-1)
#             for i in lines.splitlines():
#                 DataID = i.split(";")[0].split("_")
#                 Val = i.split(";")[1]
#     #             print("DataID", DataID)
#     #             print("Val", Val)
#                 if "LeafOBJ" in DataID[0]:
#                     SlotVal = int(DataID[1])*NumberOfLeafes + int(DataID[2])
#                     ValLeafList[SlotVal] = Val #One value per segment instead of one value per quadrilateral like in Irradiance_Leaf
#     #                 print("SlotVal", SlotVal)
#                     for l in range(QuadIdList[SlotVal], QuadIdList[SlotVal+1]):
#     #                     print("SlotVal", SlotVal,"Number of Quadrilaterals",vp.QuadCounterListList[SlotVal],"Range",QuadIdList[SlotVal],QuadIdList[SlotVal+1],"Actual writting", l)
#                         Irradiance_Leaf[l] = float(Val)*500000 # !!!!!! REMOVE 500000 !!!!! for testing purpose.
#                 if "SegsOBJ" in DataID[0]:
#                     Irradiance_Stem[DataID[1]] = float(Val) #one value per segment

#         f.close()
        ValLeafGlobal = [0]*SegsLength
#         for k, i in enumerate(LeafSegIDList):
#             ValLeafGlobal[i] = ValLeafList[k]
    #     print("ValLeafList",len(ValLeafList),ValLeafList, "\n") #data in organ IDs order (can have multiple times the same ID between organs)
    #     print("ValLeafGlobal",len(ValLeafGlobal),ValLeafGlobal, "\n") #data in organism IDs order or global IDs. (each ID in the organism is unique)
    #     print("Irradiance_Stem", len(Irradiance_Stem),Irradiance_Stem, "\n") #values already linked to the global IDs


        Multiplier = 1
        NormalsZValue = [element * Multiplier for element in NormalsZValue]
        AngleVector = [0]*len(NormalsZValue)
        for c,element in enumerate(NormalsZValue):
            AngleVector[c] = math.degrees(math.acos(element))
#         print("AngleVector", len(AngleVector), AngleVector)
        
        NormalsZValueLeaf = [0]*QuadIdList[-1]
        AngleVectorLeaf = [0]*QuadIdList[-1]
        for i in range(len(QuadIdList)-1):
            for l in range(QuadIdList[i], QuadIdList[i+1]):
                NormalsZValueLeaf[l] = NormalsZValue[i]
                AngleVectorLeaf[l] = AngleVector[i]
#         print("NormalsZValue",len(NormalsZValue),NormalsZValue)
#         print("AngleVectorLeaf", len(AngleVectorLeaf), AngleVectorLeaf, "\n")
#         print("NormalsZValueLeaf",len(NormalsZValueLeaf),NormalsZValueLeaf, "\n")
        
        
        
        #Note: the next line doesn't work if the code didn't go through a "plot_plant" first to create the "vp.QuadCounterListList"
        vp.plot_plant(plant, "NormalsZValue", ExtraParam = [ ["Irradiance", Test_Stem, Test_Leaf] , ["NormalsZValue", StemNormalsZValue, NormalsZValueLeaf], ["Angle",StemNormalsZValue,AngleVectorLeaf] ])
        #example of usage of ExtraParam:
#         names = ["Test","Angle"]
#         for name in names:
#             vp.plot_plant(plant, name, ExtraParam = [ ["Irradiance", Test_Stem, Test_Leaf] , ["Test", StemNormalsZValue, NormalsZValueLeaf], ["Angle",StemNormalsZValue,AngleVectorLeaf] ])
    
    
    # zoom on leaf--theory--2D
    print("2D leaf shape of a full grown leaf")
    lorg = plant.getOrgans(pb.leaf)[0]
    lrp = lorg.getLeafRandomParameter()    
    
    if leafRadial:
        N = 100
        yy = np.linspace(0, lorg.leafLength(), N)
        geom_x, geom_y = [],[]
        for i, x in enumerate(lrp.leafGeometry):
            geom_x.extend(x)
            geom_y.extend([yy[i]] * len(x))
        geom_x = np.array(geom_x)
        geom_y = np.array(geom_y)        
        a  = lorg.leafArea() / lorg.leafLength() # scale radius
        plt.plot(geom_x * a, geom_y, "g*")
        plt.plot(-geom_x * a, geom_y, "g*")
    else:
        geom_x_a =  np.array([0])
        geom_x_b = np.array([ x[-1] for x in lrp.leafGeometry]) #normalized x value along length
        geom_x = np.concatenate((geom_x_a,geom_x_b))
        geom_y_a = np.array([0])
        geom_y_b =np.linspace(lrp.lb, lorg.leafLength()+lrp.lb, len(geom_x_b))
        geom_y = np.concatenate((geom_y_a,geom_y_b))
        a  = lorg.leafArea() / lorg.leafLength() # scale radius
        plt.plot(geom_x * a, geom_y, "g-*")
        plt.plot(-geom_x * a, geom_y, "g-*")
    plt.ylim([0, lrp.lmax+1])
    plt.xlim([-a-1, a+1])
    plt.axis('scaled')
    plt.show()
    

r = Leuning(plant) 
r.setKr([[kr],[kr_stem],[gmax]]) #gmax will be changed by the leuning function 
r.setKx([[kz]])
leaf_nodes = r.get_nodes_index(4)

# Numerical solution 
results=[]
resultsAn=[]
resultsgco2=[]
resultsVc=[]
resultsVj=[]
resultscics=[]
resultsfw=[]
resultspl=[]

# for i in range(t_init, t_end):
if True:
    i = t_init
    print(i)
    Q_input = df['PAR'][i]
#     print("Q_input",Q_input, type(Q_input))
    RH_input = df['RH'][i]
    Tair_input = df['Tair'][i]
    p_s_input = df['ps'][i]
    print((p_s_input))
    N_input = 4.4 #nitrogen satisfaction for small wheat plants
    cs_input = df['co2'][i]
    
    p_top = p_s_input
    p_bot = p_top - 35
    p_s = np.linspace(p_top, p_bot, 35)
    
    Rgaz=8.314 #J K-1 mol-1 = cm^3*MPa/K/mol
    rho_h2o = 1 #g/cm3
    Mh2o = 18.05 #g/mol
    MPa2hPa = 10000
    hPa2cm = 1.0197
    TairK = Tair_input + 273.15
    p_a = np.log(RH_input)*Rgaz*rho_h2o*TairK/Mh2o *MPa2hPa*hPa2cm #air water potential [cm]

    r.airPressure = p_a
    
    
    
    Multiplier = Q_input
    
    tempNormalsZValue = [0]*len(NormalsZValue)
    for count, element in enumerate(NormalsZValue):
        tempNormalsZValue[count] = element * Multiplier
#         print(type(element))
#     print(type(NormalsZValue))
    tempNormalsZValue = np.array(tempNormalsZValue)
#     print(tempNormalsZValue)
#     var = [tempNormalsZValue, RH_input, TairK, p_s, N_input, cs_input]
    var = [Q_input, RH_input, TairK, p_s, N_input, cs_input]
    
    
    es = 0.61078 * np.exp(17.27 * var[2] / (var[2] + 273.3)) #FAO56
    ea = es * var[1]
    VPD = es - ea 
    r.Param['Patm']=df['Pair'][i]
    rx = r.solve_leuning( sim_time = simtime,sxx=var[3], cells = True, Qlight = var[0],VPD = VPD,Tl = var[2] + 273.15,
                         p_linit = p_s_input,
    ci_init = var[5]*0.7, cs =var[5], soil_k = [], N = var[4], log = False, verbose = False)
    fluxes = r.radial_fluxes(simtime, rx, var[3], k_soil, True)  # cm3/day
    organTypes = np.array(r.rs.organTypes)  #cm³/jour
    results.append(sum(np.where(organTypes == 4, fluxes,0)))
#     resultsAn.append(np.mean(r.An)*1e6)
    for elem in r.An:
        resultsAn.append(elem*1e6)
#     resultsAn.append((r.An)*1e6)
#     resultsVc.append(np.mean(r.Vc)*1e6)
    for elem in r.Vc:
        resultsVc.append((elem)*1e6)
#     resultsVj.append(np.mean(r.Vj)*1e6)
    for elem in r.Vj:
        resultsVj.append(elem*1e6)
#     resultsgco2.append(np.mean(r.gco2))
    for elem in r.gco2:
        resultsgco2.append((elem))
#     resultscics.append(np.mean(r.ci)/var[5])
    for elem in ((r.ci)/var[5]):
        resultscics.append(elem)
#     resultsfw.append(np.mean(r.fw))
    for elem in r.fw:
        resultsfw.append(elem)
#     resultspl.append(np.mean(r.x[leaf_nodes]))
    for elem in r.x[leaf_nodes]:
        resultspl.append(elem)
    print("resultsAn", len(resultsAn),"\n","resultsVc", len(resultsVc),"\n","resultsVj", len(resultsVj),"\n","resultsgco2", len(resultsgco2),"\n","resultscics", len(resultscics),"\n","resultsfw", len(resultsfw),"\n","resultspl", len(resultspl),"\n")

# vp.plot_plant(plant, "Test", ExtraParam = [ ["Irradiance", Test_Stem, Test_Leaf] , ["Test", StemNormalsZValue, NormalsZValue] ])


print("hey")
# plot results 
nodes = r.get_nodes()
fig, ax = plt.subplots()
name = ["root", "stem", "leaf"]
color = ['tab:blue', 'tab:orange', 'tab:green']
for ndType in [2, 3, 4]:
    y = r.get_nodes_organ_type(ndType)#coordinates
    x = rx[r.get_nodes_index(ndType)]
    ax.scatter(x, y[:,2], c=color[ndType-2],  label=name[ndType-2],
               alpha=0.3, edgecolors='none')

ax.legend()
ax.grid(True)
plt.xlabel("Xylem pressure (cm)")
plt.ylabel("Depth (cm)")
plt.title("Xylem matric potential (cm)")
plt.show()


fig, ax = plt.subplots()
name = ["root", "stem", "leaf"]
color = ['tab:blue', 'tab:orange', 'tab:green']
for ndType in [2, 3, 4]:
    segIdx = r.get_segments_index(ndType)
    nodesy = segIdx + np.ones(segIdx.shape, dtype = np.int64)
    y = nodes[nodesy]#coordinates
    x = fluxes[segIdx]
    ax.scatter(x, y[:,2], c=color[ndType-2],  label=name[ndType-2],
               alpha=0.3, edgecolors='none')

ax.legend()
ax.grid(True)
plt.xlabel("Fluxes (cm3/day)")
plt.ylabel("Depth (cm)")
plt.title("water fluxes")
plt.show()