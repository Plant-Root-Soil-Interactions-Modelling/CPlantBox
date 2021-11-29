import sys; sys.path.append("../../src/python_modules/"); sys.path.append("../../")

import plantbox as pb
import rsml_reader
from rsml_data import RsmlData
import xylem_flux

import numpy as np


class EstimateDataModel:
    """ 
    Model in the sense of model view controller (MVC), stores most that is presented in the view
    
    """

    def __init__(self):
        rsml_data = []  # list of rsml data to analyse, rsml_data.RsmlData

    def open_folder(self, folder_name):
        """ see RsmlData.open_rsml() in src/python_modules/rsml_data.py                
        """
        pass  # TODO

