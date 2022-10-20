""" water movement within the root (static soil) """


import sys; 
#from joblib import Parallel, delayed
import math
import os
import numpy as np

from datetime import datetime, timedelta
#from husain_i import runSim

main_dir=os.environ['PWD']#dir of the file
directoryN = "husain_i_TEST"#"/"+os.path.basename(__file__)[:-3]+"/"#"/a_threshold/"
results_dir = main_dir +"/results"+directoryN


if not os.path.exists(results_dir):
    os.makedirs(results_dir)
else:
   test = os.listdir(results_dir)
   for item in test:
        try:
            os.remove(results_dir+item)
        except:
            pass


