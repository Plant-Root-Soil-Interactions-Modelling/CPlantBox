""" map root segments to a soil grid """
import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import visualisation.vtk_plot as vp

import numpy as np

""" parameters """
sim_time = 14  # [day]
rs_age = 1  # initial age
age_dependent = False  # conductivities
dt = 0.1  # [days] Time step must be very small
periodic = False

""" root system """
rs = pb.MappedRootSystem()
path = "../../modelparameter/structural/rootsystem/"
name = "Anagallis_femina_Leitner_2010"  # Zea_mays_1_Leitner_2010
rs.readParameters(path + name + ".xml")

""" soil """
min_ = np.array([-5, -5, -15])
max_ = np.array([5, 5, 0.])
res_ = np.array([1, 3, 5])
if not periodic:
    sdf = pb.SDF_PlantBox(0.99 * (max_[0] - min_[0]), 0.99 * (max_[1] - min_[1]), 0.99 * (max_[2] - min_[2]))
    rs.setGeometry(sdf)
rs.setRectangularGrid(pb.Vector3d(min_), pb.Vector3d(max_), pb.Vector3d(res_), False)  # cut and map segments

rs.initialize()
rs.simulate(rs_age, False)
N = round(sim_time / dt)

anim = vp.AnimateRoots(rs)
anim.min = min_
anim.max = max_
anim.res = res_
anim.file = "results/example7a"
anim.start()

for i in range(0, N):

    rs.simulate(dt, False)

    """ add segment indices """
    segs = rs.segments
    x = np.zeros(len(segs))
    for i, s in enumerate(segs):
        try:
            x[i] = rs.seg2cell[i]
        except:  # in case the segment is not within the domain
            x[i] = -10

    ana = pb.SegmentAnalyser(rs.mappedSegments())
    ana.addData("linear_index", x)

    anim.rootsystem = ana
    anim.root_name = "linear_index"
    anim.update()
