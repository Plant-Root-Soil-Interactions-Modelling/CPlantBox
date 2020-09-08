""" Shehan's soil core example """

import sys;  sys.path.append("../..")
import plantbox as pb

print("hallo")
import math
import numpy as np

path = "../../modelparameter/rootsystem/"
name = "wheat"  # Zea_mays_4_Leitner_2014"


def initialize_root_systems(N :int, M :int, distN :float, distM :float):
    """ 
    Initializes M*N root systems 
        
    @param N         number of rows
    @param M         number of columns
    @param distN     distance between rows
    @param distM     distance between columns  
    @return a list of initialized root systems
    """
    allRS = []
    for i in range(0, N):
        for j in range(0, M):
            rs = pb.RootSystem()
            rs.readParameters(path + name + ".xml")
            rs.getRootSystemParameter().seedPos = pb.Vector3d(distN * i, distM * j, -3.)  # cm
            rs.initialize(False)  # verbose = False
            allRS.append(rs)
    return allRS


def simulate_rs(times :list, allRS):
    """ 
    Simulates all root systems for 

    @param times     simulation times (days)
    @param allRS     list of root systems to simulate
    """
    t = times.copy()
    t.insert(0, 0.)
    dt_ = np.diff(np.array(t))
    for dt in dt_:
        for rs in allRS:
            rs.simulate(dt)


def get_result(allRS :list, time :float):
    """ 
    Retrieves a the state of the root systems at a certain time 
    in a single SegmentAnalyser object
    
    @param allRS     list of root systems 
    @param time      of the simulation result (days)
    @return SegmentAnalyser object conaining the segments of 
    all root sytstems at the specific time 
    """
    a = pb.SegmentAnalyser()
    for rs in allRS:
        a_ = pb.SegmentAnalyser(rs)
        a_.filter("creationTime", 0., time)
        a_.pack()
        a.addSegments(a_)
    return a


def soil_cores(x :list, y :list, r :float, h :float):
    """
    A lsit of soil core geometries with a fixed location in the field  
 
    @param x     x coordinates of the soil cores (cm)
    @param y     y coordinates of the soil cores (cm)
    @param r     radius of the soil core (cm)
    @param h     height of the soil core (cm)
    """
    assert len(x) == len(y), "coordinate length must be equal"
    core = pb.SDF_PlantContainer(r, r, h, False)
    cores = []
    for i in range(0, len(x)):
        cores.append(pb.SDF_RotateTranslate(core, 0., pb.SDF_Axis.xaxis, pb.Vector3d(x[i], y[i], 0.)))  # just translate
    return cores;


allRS = initialize_root_systems(3, 3, 17, 45)

times = [30, 60, 90]
simulate_rs(times, allRS)

x = [ 0, 1, 2 ]
y = [ 0, 1, 2 ]
r = 2.1  # core radius (cm)
h = 160.  # core length (cm)
cores = soil_cores(x, y, r, h)

exportVTP = True  # export single root system cores (for debugging)

dz = 5  # layer thickness (cm)
nol = round(h / dz)  # number of layers
result_matrix = np.zeros((len(times) * len(cores), nol))

for i, t in enumerate(times):

    print("Analyse time", t)

    for j in range(0, len(cores)):

        core_analyser = get_result(allRS, t)
        core_analyser.write("wheat.vtp");
        core_analyser.crop(cores[i]);
        core_analyser.pack()

        tl = core_analyser.distribution("length", 0, -h, nol, True)  # vertical length distribution
        tl = np.array(tl) / (len(cores) * r * r * math.pi * dz)  # <<< TODO WHY len(cores), its croped to 1 cylinder

        result_matrix[i * len(cores) + j, :] = tl

        if exportVTP:
            vtp_name = name + "_core_cropped" + str(i) + ".vtp";  # export cropped segments for vizualisaten
            core_analyser.write(vtp_name);

np.savetxt(name + "_core_matrix.txt", result_matrix, delimiter = ', ')

if exportVTP:
    g_name = name + "_core.py";
    cores_union = pb.SDF_Union(cores)
    allRS[0].setGeometry(cores_union)  # just for writing
    allRS[0].write(g_name)

print("fin.")
