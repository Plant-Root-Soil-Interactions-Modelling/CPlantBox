"""small example"""
import sys
sys.path.append("../../.."); sys.path.append("../../../src/python_modules")
import plantbox as pb
import vtk_plot as vp

pl = pb.MappedPlant() #pb.MappedRootSystem() #pb.MappedPlant()

path = "../../../modelparameter/plant/" #"../../../modelparameter/rootsystem/" 
name = "smallPlant_mgiraud"#"manyleaves"#"oneroot_mgiraud" #"manyleaves"
pl.readParameters(path + name + ".xml")
# Initialize
pl.initialize()
start = 1
# Simulate


pl.simulate(start, True)
phl = PhloemFluxPython(pl)
segs = pl.getPolylines() #segments regrouped per organ
phl.mp2mesh(segs) #creates grid
organTypes = phl.get_organ_types()
print(phl.mesh.length, "\n", phl.orgID, "\n", phl.orgLength, "\n", [organTypes[phl.new2oldNodeID[xi] - 1] for xi in phl.mesh.cellFaceIDs[1]])

print("2nd sim")
pl.simulate(start, True)
print("\n2nd sim end")
phl = PhloemFluxPython(pl)
phlbis = PhloemFluxPythonbis(pl)
segs = pl.getPolylines() #segments regrouped per organ
phl.mp2mesh(segs) #creates grid
organTypes = phl.get_organ_types()
print(phl.mesh.length, "\n", phl.orgID, "\n", phl.orgLength, "\n", [organTypes[phl.new2oldNodeID[xi] - 1] for xi in phl.mesh.cellFaceIDs[1]])


print("3rd sim")
pl.simulate(start, True)
print("\n3rd sim end")
phl = PhloemFluxPython(pl)
phlbis = PhloemFluxPythonbis(pl)
segs = pl.getPolylines() #segments regrouped per organ
phl.mp2mesh(segs) #creates grid
organTypes = phl.get_organ_types()
print(phl.mesh.length, "\n", phl.orgID, "\n", phl.orgLength, "\n", [organTypes[phl.new2oldNodeID[xi] - 1] for xi in phl.mesh.cellFaceIDs[1]])
print(phl.mesh.cellFaceIDs)


# Plot, using vtk
#vp.plot_roots(rs, "creationTime")
