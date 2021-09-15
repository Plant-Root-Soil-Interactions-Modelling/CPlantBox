import os
import glob

for aa in range(1,3):
    for bb in range(1,3):
        for cc in range(1,22):
            fileList = glob.glob('simulations/plant'+str(aa)+'/exudate'+str(bb)+'/day'+str(cc)+'/exud*.*', recursive=True)
            for filePath in fileList:
                try:
                    os.remove(filePath)
                except OSError:
                    print("Error while deleting file")
