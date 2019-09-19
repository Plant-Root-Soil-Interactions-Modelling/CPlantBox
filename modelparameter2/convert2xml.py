import os
from os import walk
import py_rootbox as rb

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
            rs = rb.RootSystem()
            rs.openFile(f[:-7], "")
            rs.getOrganRandomParameter(rb.OrganTypes.seed)[0].name = f[:-7]
            rs.writeParameters(f[:-7] + ".xml", "plant", comments)
            print("done.\n")

