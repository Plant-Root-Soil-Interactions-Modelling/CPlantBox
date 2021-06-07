"""small example"""
import sys
sys.path.append("../../.."); sys.path.append("../../../src/python_modules")
import plantbox as pb
import vtk_plot as vp
import numpy as np
from scipy.linalg import norm


# Open plant and root parameter from a file
pl = pb.MappedPlant() #pb.MappedRootSystem() #pb.MappedPlant()
path = "../../../modelparameter/plant/" #"../../../modelparameter/rootsystem/" 
name = "smallPlant_mgiraud"#"Anagallis_femina_Leitner_2010"
pl.readParameters(path + name + ".xml")

# Initialize
pl.initialize()

# Simulate
pl.simulate(30, True)

#post processing
segs = pl.getPolylines() 
for segnum, seg in enumerate(segs):#go through each organ
    seg = np.array(list(map(lambda x: np.array(x), seg))) # converts the list of Vector3d to a numpy array. shape: (num_nodes, 3)
    d = np.diff(seg, axis=0)
    length_segments = np.array([round(norm(di),6) for di in d])
    print(seg)
    print(length_segments)
            

# Export final result (as vtp)
pl.write("results/example_1a.vtp")

# Plot, using vtk
vp.plot_roots(pl, "creationTime")
