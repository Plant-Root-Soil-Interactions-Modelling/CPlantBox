import timeit


import numpy as np
#import matplotlib.pyplot as plt

import plantbox as pb
from plantbox import Photosynthesis
from functional.xylem_flux import XylemFluxPython 

#import rsml_reader as rsml



class PhotosynthesisPython(Photosynthesis, XylemFluxPython):
    """  wrapper for photosynthesis
       
    """

    def __init__(self, plant_, psiXylInit, ciInit):
        """ @param rs is either a pb.MappedRootSystem, pb.MappedSegments, or a string containing a rsml filename"""
        Photosynthesis.__init__(self, plant_, psiXylInit, ciInit)
        XylemFluxPython.__init__(self,plant_)

    