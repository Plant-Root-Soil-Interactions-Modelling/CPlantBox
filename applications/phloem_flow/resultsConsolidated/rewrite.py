import sys; 
import os
import numpy as np
import time
import shutil



def doRewrite(srcdir_old, results_dir):
    if not os.path.exists(srcdir_old):
        raise Exception

    if not os.path.exists(results_dir):
        os.makedirs(results_dir)
    else:
        shutil.rmtree(results_dir)
        os.makedirs(results_dir)

    allfiles = os.listdir(srcdir_old)

    for log in allfiles:
        if "input" in log:
            shutil.copy2(srcdir_old + "/"+log, results_dir + "/"+log)
            
            # with open(srcdir_old + "/"+log, 'r') as f:
            #     lines = f.read().splitlines()
            #     #last_line = lines[-1]
            # with open(results_dir + "/"+log, 'w') as f:
            #     f.write(lines[0])
            #     f.write(lines[1])
        else:
            if log.endswith(".csv"):
                with open(srcdir_old + "/"+log, 'r') as f:
                    lines = f.read().splitlines()
                    last_line = lines[-1]
                with open(results_dir + "/"+log, 'w') as f:
                    f.write(last_line)


srcdir_olds = "CalibPart2_"
for i in range(3):
    print("../results/"+srcdir_olds + repr(i +1 ), "../results/"+srcdir_olds + repr(i +1) + "_rewrote")
    doRewrite("../results/"+srcdir_olds + repr(i +1 ), "../results/"+srcdir_olds + repr(i +1) + "_rewrote" )
    