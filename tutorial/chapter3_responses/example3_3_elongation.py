"""scales the root elongation rate"""

import plantbox as pb
import plantbox.visualisation.vtk_plot as vp

plant = pb.Plant()
path = "../../modelparameter/structural/rootsystem/"
filename = "Anagallis_femina_Leitner_2010"
plant.readParameters(path + filename + ".xml")

# box with a left and a right compartment for analysis
sideBox = pb.SDF_PlantBox(10, 20, 50)
left = pb.SDF_RotateTranslate(sideBox, pb.Vector3d(-4.99, 0, 0))
right = pb.SDF_RotateTranslate(sideBox, pb.Vector3d(4.99, 0, 0))
leftright = pb.SDF_Union(left, right)
plant.setGeometry(leftright)

# left compartment has a minimum of 0.01, 1 elsewhere
maxS = 1.0  # maximal
minS = 0.05  # minimal
slope = 1.0  # [cm] linear gradient between min and max
leftC = pb.SDF_Complement(left)
soilprop = pb.SoilLookUpSDF(leftC, maxS, minS, slope)

# Manually set scaling function and tropism parameters
for organ_type in [pb.root, pb.stem, pb.stem]:
    sigma = [0.4, 1.0, 1.0, 1.0, 1.0] * 2
    for p in plant.getOrganRandomParameter(organ_type):
        p.dx = 0.25  # adjust resolutionx
        p.tropismS = sigma[p.subType - 1]
        p.f_se = soilprop  # 1. Scale elongation

# simulation
plant.initialize()
sim_time = 60.0
dt = 1.0
for i in range(0, round(sim_time / dt)):
    # in a dynamic setting change soilprop here
    plant.simulate(dt, False)

# analysis
l = plant.getSummed("length")
al = pb.SegmentAnalyser(plant)
al.crop(left)
ll = al.getSummed("length")
ar = pb.SegmentAnalyser(plant)
ar.crop(right)
lr = ar.getSummed("length")
print(f"\nLeft  compartment total root length {ll:g} cm, {100 * ll / l:g}%")
print(f"\nRight compartment total root length {lr:g} cm, {100 * lr / l:g}% \n")

# write results
plant.write("results/example_5a.py")  # compartment geometry
plant.write("results/example_5a.vtp")  # root system

# plot, using vtk
vp.plot_roots(plant, "rootLength")  # press 'y'
