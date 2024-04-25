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
    adds functionality to initialize with a static root system, where emerging lateral roots are predefined
    static root system has subType (57 for ST)
    
    usage:    
    PlantPython.initialize_static (instead of initialize), no organs emerge from the plant seed
    define lateral types and emergence times
    PlantPython.initialize_static_laterals() to create the dynamic organs    
    """

    def __init__(self, seed_num = 0):
        super().__init__(seed_num)
        self.static_organs = {}
        self.data = None  # RsmlData

    def plot_rsml_(self, polylines:list, prop:list):
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

    def polyline_length_(self, i0, i1, poly):
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
        radii, et, types, tag_names, _, _ = rsml_reader.get_root_parameters(self.data.polylines, self.data.functions, self.data.properties)
        # types = self.data.properties[self.data.tagnames[2]]  # e.g 'order'
        # print(self.data.tagnames)
        # print(self.data.properties.keys(), len(self.data.polylines))

        seed = Seed(self)

        # 2. Parse rsml, create self.static_organs
        parent_id = {}
        static_polys, static_types = [], []
        for i, root in enumerate(self.data.polylines):
            if types[i] in initial_sub_types:
                parent_id = self.data.properties["parent-poly"][i]
                pni = self.data.properties["parent-node"][i]

                length = self.polyline_length_(0, len(root) - 1, root)
                a = np.mean(radii[i])
                # print(a, len(radii[i]))
                r = 1.e6  # planted initially, and static
                lb = 0  # length of basal zone
                theta = 0  # insertion angle
                rlt = 1.e6  # root life time
                ln_ = []
                laterals = False
                param = RootSpecificParameter(0, lb, length, ln_, r, a, theta, rlt, laterals)  ############## which subType
                # print(param)

                id = self.getOrganIndex()  # next index
                organ = StaticRoot(id, param, length, pni)

                organ.setOrganism(self)  # needed for adding nodes
                self.static_organs[i] = organ

                static_polys.append(self.data.polylines[i])  # for debugging
                static_types.append(types[i])

        # self.plot_rsml_(static_polys, static_types)

        # 3. Add geometry
        for i, root in enumerate(self.data.polylines):
            # print("root", i, ":", len(self.data.polylines[i]), self.data.polylines[i])
            if types[i] in initial_sub_types:
                parent_id = self.data.properties["parent-poly"][i]
                pni = self.data.properties["parent-node"][i]
                organ = self.static_organs[i]
                if parent_id >= 0:
                    # print(i, "parent", parent_id, "pni", pni)
                    parent = self.static_organs[parent_id]
                    organ.addNode(parent.getNode(pni), parent.getNodeId(pni), 0.)
                    # print(self.data.polylines[parent_id][pni])
                for node in self.data.polylines[i]:
                    organ.addNode(Vector3d(node), 0.)

        # 4. Create topology of static roots
        for i, root in enumerate(self.data.polylines):
            if types[i] in initial_sub_types:
                parent_id = self.data.properties["parent-poly"][i]
                organ = self.static_organs[i]
                try:
                    parent = self.static_organs[parent_id]
                    # organ.setParent(parent)  # actually in addChild
                    parent.addChild(organ)
                    # print("Added", i, "to parent", parent_id)
                except:
                    print("PlantPython: initialize_static(): organ", i, "has no parent", parent_id)
                    seed.addChild(self.static_organs[0])
        self.addOrgan(seed)

        # 4. The CPlantBox part
        # seed.initialize() # not called i.e. no tap root or basal roots are created
        self.oldNumberOfNodes = self.getNumberOfNodes()
        self.initCallbacks()

    def set_identical_laterals(self, initial_sub_types, lateral_subtypes, emerge_type):
        """ places laterals as in the original rsml, all start growing at once 
        
        initial_sub_types             subTypes of the initial static root system 
        lateral_subtypes              subTypes of the laterals wihtin the RSML
        emerge_type                   subType of the model lateral (todo: could be probabilistic)        
        """
        lni, lt, ld = [], [], []
        types = self.data.properties[self.data.tagnames[2]]
        for i, root in enumerate(self.data.polylines):
            if types[i] in lateral_subtypes:
                parent_id = self.data.properties["parent-poly"][i]
                if types[parent_id] in initial_sub_types:
                    pni = self.data.properties["parent-node"][i] + 1
                    parent = self.static_organs[parent_id]
                    parent.addLateral(pni, emerge_type, 0.)
                    parent.param().laterals = True
                    # print("Root", i, ":", parent_id, pni, emerge_type, 0.)

    def initialize_static_laterals(self):
        """ Creates lateral root instances from static roots lateralDelays, lateralTypes, and lateralDelays """
        for organ in self.static_organs.values():
            organ.initializeLaterals()

    def analyse_laterals(self, initial_sub_types, lateral_subtypes):
        """ prints emergence points of laterals (for debugging)"""
        print("Plant.analyse_laterals()")
        laterals, tip_laterals = {}, {}
        types = self.data.properties[self.data.tagnames[2]]
        for i, root in enumerate(self.data.polylines):
            if types[i] in lateral_subtypes:
                parent_id = self.data.properties["parent-poly"][i]
                if types[parent_id] in initial_sub_types:
                    pni = self.data.properties["parent-node"][i]
                    if parent_id in laterals:
                        laterals[parent_id].append(i)
                    else:
                        laterals[parent_id] = [i]
                    tip_lateral = pni == len(self.data.polylines[parent_id]) - 1
                    if tip_lateral:
                        if parent_id in tip_laterals:
                            tip_laterals[parent_id].append(i)
                        else:
                            tip_laterals[parent_id] = [i]
                    print("Id", i, "subType", types[i], "parent", parent_id, "parent subType", types[parent_id], "at node index", pni, "located at tip", tip_lateral)

        return laterals, tip_laterals


if __name__ == '__main__':

    plant = PlantPython()

    # path = "../../modelparameter/structural/rootsystem/"
    # name = "Glycine_max_Moraes2020"
    # plant.readParameters(path + name + ".xml")
    plant.readParameters("wine.xml")
    # plant.setOrganRandomParameter(SeedRandomParameter(plant))

    p2 = plant.getOrganRandomParameter(2, 2)
    p2.successor = [[3]]
    p2.successorP = [[1]]
    print(p2)
    # p2.theta = 0
    # p2.thetas = 0
    p2.tropismT = 1  #  1 gravi, 2 exo
    p2.tropismS = 0.2
    p2.tropismN = 0.5  # 0.05

    plant.initialize_static("B-23.rsml", [0, 1])  # 0 is shoot, 1 are static roots

    print()
    laterals, tip_laterals = plant.analyse_laterals([0, 1], [2, 3])  # for debugging
    print()

    plant.set_identical_laterals([0, 1], [2, 3], 2)
    plant.initialize_static_laterals()

    plant.simulate(125., True)

    vp.plot_roots(plant, "subType")

    print("fin")
