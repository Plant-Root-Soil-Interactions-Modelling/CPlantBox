"""Scale root elongation based on EquidistantGrid1D"""

import numpy as np
import plantbox as pb
import plantbox.visualisation.vtk_plot as vp

rs = pb.Plant()
path = "../../modelparameter/structural/rootsystem/"
name = "Anagallis_femina_Leitner_2010"
rs.readParameters(path + name + ".xml")

scale_elongation = pb.EquidistantGrid1D(0, -50, 100)  # |\label{l3_4_scale:gridStart}|
soil_strength = np.ones((99,))
soil_strength[20:30] = 5  # |\label{l3_4_scale:gridEnd}|

scales = np.exp(-0.4 * soil_strength)  #  exemplary scaling function # |\label{l3_4_scale:functionStart}|
scale_elongation.data = scales  # set proportionality factors
print("-5 cm ", scale_elongation.getValue(pb.Vector3d(0, 0, -5)))
print("-15 cm", scale_elongation.getValue(pb.Vector3d(0, 0, -15)))  # |\label{l3_4_scale:functionEnd}|

for organ_type in [pb.root, pb.stem, pb.leaf]:
    for p in rs.getOrganRandomParameter(organ_type):
        p.f_se = scale_elongation  # |\label{l3_4_scale:applyScaling}|

rs.initialize()
simtime = 30.
dt = 0.1  # small, for animation

anim = vp.AnimateRoots(rs)
anim.root_name = "creationTime"
anim.file = "example3_4_scale_elongation"

# outlines dense layer in white
anim.min = np.array([-10, -10, -15])
anim.max = np.array([10, 10, -10.])
anim.res = np.array([1, 1, 1])
anim.start()

for i in range(0, round(simtime / dt)):  # Simulation

    # option for dynamic scale update
    new_scales = scales  # |\label{l3_4_scale:dynamicStart}|
    scale_elongation.data = new_scales  # |\label{l3_4_scale:dynamicEnd}|

    rs.simulate(dt, False)
    anim.update()

rs.write("results/example3_4_scale_elongation.vtp")
anim.iren.Start()  # Keeps the render window open
