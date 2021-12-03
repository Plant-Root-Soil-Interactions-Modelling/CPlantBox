import sys; sys.path.append("../../src/python_modules/"); sys.path.append("../../")

from rsml_data import RsmlData
import plantbox as pb
import rsml_reader

import numpy as np
import os


class EstimateDataModel:
    """ 
    Model in the sense of model view controller (MVC), stores most that is presented in the view    
    """

    def __init__(self):
        self.rsmls = []  # list of rsml data to analyse, rsml_data.RsmlData
        self.times = []
        self.base_root_indices = []
        self.folder_name = ""

    def open_folder(self, folder_name):
        """ see RsmlData.open_rsml() in src/python_modules/rsml_data.py                
        """
        self.rsmls = []  # delete old data
        self.folder_name = folder_name
        for root, subdirs, files in os.walk(folder_name):
            for filename in files:
                if filename.endswith('.rsml'):
                    file_path = os.path.join(root, filename)
                    print('file %s (full path: %s)\n' % (filename, file_path))
                    file_data = RsmlData()
                    file_data.open_rsml(file_path)
                    self.rsmls.append(file_data)
                    try:
                        str = ''.join([n for n in filename if n.isdigit()])
                        self.times.append(int(str))
                    except:
                        self.times.append(0)
        print("times", self.times)
        self.base_roots_()  # find base roots
        self.create_length()  # add a length tag to properties

    def exists(self):
        """ true if a rsml folder was set """
        return len(self.rsmls) > 0

    def base_roots_(self):
        """ return only the base roots """
        self.base_root_indices = []
        for i, data in enumerate(self.rsmls):
            pp = data.properties["parent-poly"]
            self.base_root_indices.append([])
            for j, _ in enumerate(data.polylines):
                if pp[j] == -1:  # add if base root
                    self.base_root_indices[-1].append(j)
        print(self.base_root_indices)

    def create_length(self):
        """ calculates root lengths, and adds it to properties if not already present)"""
        for k, _ in enumerate(self.rsmls):
            if not "length" in self.rsmls[k].properties:
                self.rsmls[k].properties["length"] = []
                for i, p in enumerate(self.rsmls[k].polylines):
                    l = 0.
                    for j, n1 in enumerate(p[:-1]):
                        n2 = p[j + 1]
                        l += np.linalg.norm(np.array(n2) - np.array(n1))
                    self.rsmls[k].properties["length"].append(l)
        else:
            print("EstimateDataModel.create_length: 'length' tag is already available")

