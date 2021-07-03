"""small example"""
import sys
sys.path.append("../../.."); sys.path.append("../../../src/python_modules")
import plantbox as pb
import vtk_plot as vp
import numpy as np
from scipy.linalg import norm


# Open plant and root parameter from a file
pl = pb.Plant() #pb.MappedRootSystem() #pb.MappedPlant()
path = "../../../modelparameter/plant/" #"../../../modelparameter/rootsystem/" 
name =  "Triticum_aestivum_adapted_2021"#"Heliantus_Pagès_2013"#"morning_glory_7m"#"Anagallis_femina_Leitner_2010"##"Anagallis_femina_Leitner_2010"#"smallPlant_mgiraud"#
pl.readParameters(path + name + ".xml")

# Initialize
pl.initialize()
# Simulate
pl.simulate(15, False)


segs = pl.getPolylines() 
lengthtot = 0
for segnum, seg in enumerate(segs):#go through each organ
    seg = np.array(list(map(lambda x: np.array(x), seg))) # converts the list of Vector3d to a numpy array. shape: (num_nodes, 3)
    d = np.diff(seg, axis=0)
    length_segments = np.array([round(norm(di),6) for di in d])
    #print(seg)
    lengthtot += sum(length_segments)
print(lengthtot)

for o in pl.getOrgans():
    print(o.param())        

for o in pl.getOrgans():
    print(o.getLength(),o.organType(), o.getId() )         
# Export final result (as vtp)
pl.write("results/example_1e.vtp")
"""
# Simulate
lista = np.zeros(35)
for i in range(1,35):
    print("day n°",i)
    pl.simulate(1, False)

    #post processing

    segs = pl.getPolylines() 
    lengthtot = 0
    for segnum, seg in enumerate(segs):#go through each organ
        seg = np.array(list(map(lambda x: np.array(x), seg))) # converts the list of Vector3d to a numpy array. shape: (num_nodes, 3)
        d = np.diff(seg, axis=0)
        length_segments = np.array([round(norm(di),6) for di in d])
        #print(seg)
        lengthtot += sum(length_segments)
    print(lengthtot)
    lista[i] = lengthtot

print(','.join([num for num in map(str,lista)]))

for o in pl.getOrgans():
    if(o.organType()>2):
        print(o.param())        

for o in pl.getOrgans():
    print(o.getLength(True),o.getLength(False),o.organType(), o.getId() )         
# Export final result (as vtp)
pl.write("results/example_1e.vtp")

# Plot, using vtk
#vp.plot_roots(pl, "creationTime")
"""