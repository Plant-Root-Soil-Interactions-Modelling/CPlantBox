"""Scale root elongation based on EquidistantGrid1D"""
import sys; sys.path.append("../../..")
import plantbox as pb
import numpy as np
import vtk_plot as vp

rs = pb.RootSystem()
path = "../../../modelparameter/rootsystem/"
name = "Anagallis_femina_Leitner_2010"
rs.readParameters(path + name + ".xml")

scale_elongation = pb.EquidistantGrid1D(0, -100, 100)  # for root elongation from 0 cm to -50 cm, 100 layers
soil_strength = np.ones((99,))
soil_strength[10:20] = 5  # data, with a very dense layer
scales = np.exp(-0.4 * soil_strength)  # scales from some empirical equation (TODO)
scale_elongation.data = scales  # set proportionality factors
print("value at -3 cm ", scale_elongation.getValue(pb.Vector3d(0, 0, -3)))
print("value at -25 cm", scale_elongation.getValue(pb.Vector3d(0, 0, -25)))

for p in rs.getRootRandomParameter():  # Manually set scale elongation function
    p.f_se = scale_elongation

# Simulation
rs.initialize()
simtime = 60.
dt = 1.
for i in range(0, round(simtime / dt)):

    # update soil model (e.g. soil_strength)

    # update scales (e.g. from water content, soil_strength)
    scales = np.exp(-0.4 * soil_strength)  # scales from some equation (TODO)

    # copy scales into scaling funciton
    scale_elongation.data = scales

    rs.simulate(dt, False)

rs.write("../results/example_5b.vtp")
vp.plot_roots(rs, "age")
