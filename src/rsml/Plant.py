import sys; sys.path.append(".."); sys.path.append("../..");

import plantbox as pb
from plantbox import *

from rsml.rsml_data import RsmlData
import rsml.rsml_reader as rsml_reader
import visualisation.vtk_plot as vp

import numpy as np
import matplotlib.pyplot as plt
import matplotlib


class PlantPython(Plant):
    """
    adds functionality to initialize with a static root system
    
    emerging lateral roots are predefined
    """

    def __init__(self, seed_num = 0):
        super().__init__(seed_num)
        self.static_organs = {}
        self.data = None  # RsmlData

    def plot_rsml(self, polylines:list, prop:list):
        """Plots the polylines in y-z axis with colors given by a root property
    
        Args:
        polylines(list): flat list of polylines, one polyline per root 
        prop(list): a single property, list of scalar value, on per root 
        """
        f = matplotlib.colors.Normalize(vmin = min(prop), vmax = max(prop))
        cmap = plt.get_cmap("jet", 256)
        for i, pl in enumerate(polylines):
            nodes = np.array(pl)
            plt.plot(nodes[:, 1], nodes[:, 2], color = cmap(f(prop[i])))
        plt.axis('equal')
        plt.show()

    def polyline_length(self, i0, i1, poly):
        """ length of poly between two indices i0 and i1 """
        l = 0.
        for i in range(i0, i1):
            l += np.linalg.norm(np.array(poly[i]) - np.array(poly[i + 1]))
        return l

    def initialize_static(self, rsml_file, initial_sub_types):
        """ uses initial_sub_types of the rsml_file as a static initial root system, 
        called instead of initialize(), initializeLB(), or initializeDB() """

        self.reset()  # resets plant

        # 1. Open the RSML file
        self.data = RsmlData()
        self.data.open_rsml(rsml_file)
        # radii, cts, types, tagnames = rsml_reader.get_parameter(self.data.polylines, self.data.functions, self.data.properties)
        radii, et, types, tag_names, _, _ = rsml_reader.get_root_parameters(self.data.polylines, self.data.functions, self.data.properties)
        # types = self.data.properties[self.data.tagnames[2]]  # e.g 'order'
        # print(self.data.tagnames)
        # print(self.data.properties.keys(), len(self.data.polylines))

        # 2. Parse rsml, create self.static_organs
        parent_poly_ = {}
        static_polys, static_types = [], []
        for i, root in enumerate(self.data.polylines):
            if types[i] in initial_sub_types:
                parent_poly = self.data.properties["parent-poly"][i]
                pni = self.data.properties["parent-node"][i]

                length = self.polyline_length(0, len(root) - 1, root)

                a = np.mean(radii[i])
                print(a, len(radii[i]))

                r = 1.e6  # planted initially, and static
                lb = 0
                theta = 0
                rlt = 1.e6
                ln_ = []
                laterals = False
                param = RootSpecificParameter(types[i], lb, length, ln_, r, a, theta, rlt, laterals)
                # print(param)

                id = self.getOrganIndex()  # next index
                organ = StaticRoot(id, param, length, pni)
                organ.setOrganism(self)  # needed for adding nodes
                parent_poly_[i] = parent_poly
                self.static_organs[i] = organ

                if parent_poly >= 0:
                    organ.addNode(Vector3d(self.data.polylines[parent_poly][pni]), 0.)
                for node in self.data.polylines[i]:
                    organ.addNode(Vector3d(node), 0.)

                static_polys.append(self.data.polylines[i])  # for debugging
                static_types.append(types[i])

        # self.plot_rsml(static_polys, static_types)

        # 3. Create topology
        for i, root in enumerate(self.data.polylines):
            if types[i] in initial_sub_types:
                organ = self.static_organs[i]
                try:
                    parent = self.static_organs[parent_poly_[i]]
                    organ.setParent(parent)
                except:
                    print("PlantPython: initialize_static(): organ", i, "has no parent", parent_poly_[i])

        # 3. Add StaticRoot to self
        for i, root in enumerate(self.data.polylines):
            if types[i] in initial_sub_types:
                organ = self.static_organs[i]
                self.addOrgan(organ)

        # 4. The CPlantBox part
        seed = Seed(self)
        self.addOrgan(seed);
        # seed.initialize() # not called i.e. no tap root or basal roots are created
        self.oldNumberOfNodes = self.getNumberOfNodes()
        self.initCallbacks()

    def initialize_static_laterals(self):
        """ Creates lateral root instances from static roots lateralDelays, lateralTypes, and lateralDelays """
        for organ in self.static_organs.values():
            organ.initializeLaterals()

    def analyse_laterals(self, initial_sub_types, lateral_subtypes):
        """ prints emergence points of laterals """
        print("Plant.analyse_laterals()")
        types = self.data.properties[self.data.tagnames[2]]
        for i, root in enumerate(self.data.polylines):
            if types[i] in lateral_subtypes:
                parent_poly = self.data.properties["parent-poly"][i]
                if types[parent_poly] in initial_sub_types:
                    pni = self.data.properties["parent-node"][i]
                    print("Id", i, "subType", types[i], "parent", parent_poly, "parent subType", types[parent_poly], "at node index", pni)
                    if pni == len(self.data.polylines[parent_poly]) - 1:
                        print("tip")


if __name__ == '__main__':

    plant = PlantPython()

    plant.setOrganRandomParameter(SeedRandomParameter(plant))

    plant.initialize_static("B-23.rsml", [0, 1])  # 0 is shoot, 1 are static roots

    # print()
    # plant.analyse_laterals([0, 1], [2, 3])

    # place clever lateral emergence logic here

    plant.initialize_static_laterals()
    plant.simulate(2, True)

    vp.plot_roots(plant, "radius")

    print("fin")
