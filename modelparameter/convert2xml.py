""" converts old parameter files (.pparam, .rparam) to the new xml format """

import os
from os import walk
import sys
sys.path.append("../..")
import plantbox as pb

#
comments = True  # why not

# go to the right place
path = os.path.dirname(os.path.realpath(__file__))
os.chdir(path)

f = []
for (dirpath, dirnames, filenames) in walk(path):
    for f in filenames:
        if f[-6:] == "pparam":
            print(f[:-7])
            rs = pb.RootSystem()
            rs.openFile(f[:-7], "")
            rs.getOrganRandomParameter(pb.OrganTypes.seed)[0].name = f[:-7]
            rs.writeParameters(f[:-7] + ".xml", "plant", comments)
            print("done.\n")

