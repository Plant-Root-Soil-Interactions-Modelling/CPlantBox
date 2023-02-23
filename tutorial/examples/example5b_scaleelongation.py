"""Scale root elongation based on EquidistantGrid1D"""
import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import visualisation.vtk_plot as vp

import numpy as np

rs = pb.RootSystem()
path = "../../modelparameter/structural/rootsystem/"
name = "Anagallis_femina_Leitner_2010"
rs.readParameters(path + name + ".xml")

scale_elongation = pb.EquidistantGrid1D(0, -50, 100)
soil_strength = np.ones((99,))
soil_strength[20:30] = 5  # data, with a very dense layer at -10 to -15 cm
scales = np.exp(-0.4 * soil_strength)  #  equation (TODO)
scale_elongation.data = scales  # set proportionality factors
print("-3 cm ", scale_elongation.getValue(pb.Vector3d(0, 0, -3)))
print("-25 cm", scale_elongation.getValue(pb.Vector3d(0, 0, -25)))

for p in rs.getRootRandomParameter():
    p.f_se = scale_elongation  # set scale elongation function

rs.initialize()

anim = vp.AnimateRoots(rs)
anim.root_name = "creationTime"
anim.file = "example5b"
anim.min = np.array([-10, -10, -50])
anim.max = np.array([10, 10, 0.])
anim.res = np.array([1, 1, 1])
anim.start()
path = "../../../modelparameter/rootsystem/"
simtime = 60.
dt = 0.1  # small, for animation
for i in range(0, round(simtime / dt)):  # Simulation

    # update soil model (e.g. soil_strength)

    # update scales (e.g. from water content, soil_strength)
    scales = np.exp(-0.4 * soil_strength)  # (TODO)

    # copy scales into scaling funciton
    scale_elongation.data = scales

    rs.simulate(dt, False)

    anim.update()

rs.write("../results/example_5b.vtp")
vp.plot_roots(rs, "age")
