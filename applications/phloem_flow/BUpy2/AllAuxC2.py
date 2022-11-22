import sys; 
import os
from AllAuxCmaster import AllAuxCmasterFunc


directoryN = "/"+os.path.basename(__file__)[:-3]+"/"
AllAuxCmasterFunc(directoryN)