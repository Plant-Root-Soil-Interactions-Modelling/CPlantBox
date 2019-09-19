import sys
sys.path.append("..")
import plantbox as pb
from rb_tools import *
import numpy as np


def test_copy():

    # Parameter
    simtime = 30.  # days
    dt = 1.
    N = round(simtime / dt)  # steps
    name = "Anagallis_femina_Leitner_2010"

    print("\n1 original")
    rs = pb.RootSystem()
    rs.openFile(name)
    rs.setSeed(1)  # before initialize to mimic
    rs.initialize()
    rs2 = pb.RootSystem(rs)
    rs.simulate(simtime)
    l1 = np.sum(v2a(rs.getParameter("length")))
    print("total length", l1)

    print("\n2 copy")
    rs2.simulate(simtime)
    l2 = np.sum(v2a(rs2.getParameter("length")))
    print("total length", l2)

    print("\n3 rebuild same")
    rs3 = pb.RootSystem()
    rs3.openFile(name)
    rs3.setSeed(1)
    rs3.initialize()
    rs3.simulate(simtime)
    l3 = np.sum(v2a(rs3.getParameter("length")))
    print("total length", l3)

    return


if __name__ == '__main__':
    test_copy()
