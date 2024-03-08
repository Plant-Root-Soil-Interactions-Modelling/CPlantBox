import sys; sys.path.append(".."); sys.path.append("../..");

from rsml.rsml_data import RsmlData

import plantbox as pb
from plantbox import *

import numpy as np


class PlantPython(Plant):
    """
    adds functionality to initialize with a static root system
    
    emerging lateral roots are predefined
    """

    def polyline_length(self, i0, i1, poly):
        """ length between two indices i0 and i1 """
        l = 0.
        for i in range(i0, i1):
            l += np.linalg.norm(np.array(poly[i]) - np.array(poly[i + 1]))
        return l

    def initialize_static(self, rsml_file, initial_sub_type):

        # 1. open the file ...
        file_data = RsmlData()
        file_data.open_rsml(rsml_file)
        print(file_data.properties.keys(), len(file_data.polylines))

        # 2. parse rsml
        organ_dict = {}
        parent_poly_ = { }
        for i, root in enumerate(file_data.polylines):
            if file_data.types[i] == initial_sub_type:
                parent_poly = file_data.properties["parent-poly"][i]
                pni = file_data.properties["parent-node"][i]
                self.id = self.getOrganIndex()  # next index

                subType = file_data.types[i]
                length = self.polyline_length(0, len(root), root)
                a = np.mean(file_data.radii[i])
                r = 1.e6  # planted initially, and static
                lb = 0
                theta = 0
                rlt = 1.e6
                ln_ = []
                laterals = False
                param = RootSpecificParameter(subType, lb, length, ln_, r, a, theta, rlt, laterals)
                print(param)

                partialIHeading = 0  # don't what it is ????

                age = 0.  # because they are planted initially

                organ = StaticRoot(id, param, False, False, age, length, partialIHeading, pni, False, 0)
                parent_poly_[i] = parent_poly
                organ_dict[i] = organ

        # 3. create topology
        for i, root in enumerate(file_data.polylines):
            if file_data.types[i] == initial_sub_type:
                pass

        # 3. add StaticRoot to self


plant = PlantPython()
plant.initialize_static("B-23.rsml", 1)
