"""small example"""
import sys
sys.path.append("../../.."); sys.path.append("../../../src/python_modules")
import plantbox as pb
import vtk_plot as vp
from phloem_flux import PhloemFluxPython  

rs = pb.MappedPlant()

# Open plant and root parameter from a file
path = "../../../modelparameter/plant/"
name = "smallPlant_mgiraud"
rs.readParameters(path + name + ".xml")

# Initialize
rs.initialize()

# Simulate
rs.simulate(30, True)


# Export final result (as vtp)
rs.write("results/example_1a.vtp")
r = PhloemFluxPython(rs,psiXylInit =-500,ciInit =(350e-6)*0.5)
print(r.get_segments())
print(r.get_nodes())