""" Faba example based on measured parameters """  
import sys
sys.path.append("../../..")
import plantbox as pb
import vtk_plot as vp

import numpy as np

# a0  0.12488266186153994 cm  0.017324242340369434 cm
# a   0.05011207835516453 cm  0.006973171086743375 cm
# la  9.972678217993812 cm   1.9916938856649533 cm
# lad 5.112252737682878 days 1.0181557940876609 days
# lb  1.006696420845271 cm   0.4921345280808324 cm
# ln  0.15673771564425693 cm 0.22356340411819697 cm
# theta 116.9425337093934    46.20089933619959
# Maximal number of branches 840.0

# TODO std from fitting


class Boundary_Elongation_Impedance(pb.SoilLookUp):

    def __init__(self, geometry, impedance=0.1):
        super(Boundary_Elongation_Impedance, self).__init__()
        self.geometry = geometry
        self.impedance = impedance

    def getValue(self, pos, root):
#         print("Root", root)
#         print("Root tip position", pos)
        d = self.geometry.getDist(pos)
        if d < 0:
            return 1.
        else: 
            return 1. * self.impedance 


rs = pb.RootSystem()
p0, p1 = pb.RootRandomParameter(rs), pb.RootRandomParameter(rs)

""" tap root """
p0.name, p0.subType = "taproot", 1
p0.successor, p0.successorP = [2], [1]  # add successors
p0.dx = 0.5  # [cm] axial resolution
# from faba_taproot, faba_laterals
p0.a, p0.a_s = 0.125, 0.017  # [cm] radius TODO
p0.lb, p0.lbs = 1.007, 0.492  # [cm] basal zone 
p0.la, p0.las = 9.973, 1.992  # [cm] apical zone 
p0.lmax = 150  # [cm] maximal root length, number of lateral branching nodes = round((lmax-lb-la)/ln) + 1
p0.ln, p0.lns = 0.157, 0.224  # [cm] inter-lateral distance (16 branching nodes)  
p0.theta = 0.  # [rad]
p0.r = 2.071  # [cm/day] initial growth rate
# visual comparison
p0.tropismT, p0.tropismN, p0.tropismS = pb.TropismType.gravi, 1.8, 0.2  
 
""" laterals """
p1.name, p1.subType = "lateral", 2
p1.dx = 0.25  # [cm] axial resolution
# from faba_taproot, faba_laterals
p1.a, p1.a_s = 0.0501, 0.0069  # [cm] radius TODO 
p1.theta = (180. - 116.943) / 180.*np.pi  # [rad] 
# p1.lmax, p1.lmaxs = 3.339, 0.2 * 3.339  # fit k and r  
# p1.r, p1.rs = 3.745, 0.1 * 3.745  # fit k and r 
p1.lmax, p1.lmaxs = 50, 0.2 * 50  # fit r  
p1.r, p1.rs = 0.483, 0.1 * 0.483  # fit r
p1.r, p1.rs = 0.823, 0.1 * 0.823  # fit r (last measurement removed)
# by visual comparison
p1.tropismT, p1.tropismN, p1.tropismS = pb.TropismType.exo, 2, 0.2  

rs.setOrganRandomParameter(p0)
rs.setOrganRandomParameter(p1)

""" seed """
srp = pb.SeedRandomParameter(rs)  # with default values
srp.seedPos = pb.Vector3d(0., 0., -0.6)  # [cm] seed position TODO ?
srp.maxB = 0  # [-] number of basal roots (neglecting basal roots and shoot borne)
rs.setRootSystemParameter(srp)

""" container geometry """ 
length = 7.425  # cm
height = 14.31  # cm
rhizotron = pb.SDF_PlantBox(length, length, height)
rs.setGeometry(rhizotron)  # soilcore, or rhizotron

# elongation impedance at boundaries 
layer = 0.5  # cm
smaller_rhizotron = pb.SDF_PlantBox(length - 2 * layer, length - 2 * layer, height - layer)
scale_elongation = Boundary_Elongation_Impedance(smaller_rhizotron, 0.1)  # 1 instead of 0.1 would disable the impedance 
for p in rs.getRootRandomParameter():
    p.f_se = scale_elongation

rs.initialize()
simtimes = np.array([0, 7, 11, 15])
rs.writeParameters("faba.xml")  # for later use

T = 10  # temporal resolution 
for i, dt in enumerate(np.diff(simtimes)):
    
    for j in range(0, T * int(dt)):
        rs.simulate(1. / T, False)  #
        
    rs.write("../results/example_faba_{}_days.vtp".format(simtimes[i + 1]))

# Plot, using vtk
vp.plot_roots(rs, "creationTime")

