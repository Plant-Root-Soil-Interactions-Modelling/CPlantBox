"""
Wrapper for MappedPlant, MappedSegments, or MappedOrganism 
to add functionality that is easier implemented in Python 
"""
import sys; sys.path.append(".."); sys.path.append("../..");
sys.path.append("../rsml");
import numpy as np
from scipy import sparse
from plantbox import MappedPlant

import plantbox as pb
from plantbox import Plant
from rsml.rsml_data import RsmlData
import rsml.rsml_reader as rsml_reader

import numpy as np

class MappedPlantPython(MappedPlant):
    def __init__(self, seed_num = 0):
        super().__init__(seed_num)
        self.static_organs = {}
        self.data = None  # RsmlData
     
    def do_simulate(self, time):
        self.simulate(time)
        self.nodes_py = np.array(list(map(lambda x: np.array(x), self.getNodes())))
        self.segments_py = np.array(list(map(lambda x: np.array(x), self.getSegments())), dtype = np.int64)
            
    def toNumpy(self,cpbArray):
        """ converts the cpbArray to a numpy array """
        return np.array(list(map(lambda x: np.array(x), cpbArray)))
        
    def get_nodes(self):
        """ converts the list of Vector3d to a 2D numpy array """
        return self.nodes_py

    def get_segments(self):
        """ converts the list of Vector2i to a 2D numpy array """
        return self.segments_py

    def get_subtypes(self):
        """ segment sub types as numpy array """
        return np.array(self.subTypes)

    def get_ages(self, final_age = -1.):
        """ converts the list of nodeCT to a numpy array of segment ages
        @param final_age [day]         current root system age, (default = 0 means detect maximum from nodeCT)
        """
        cts = np.array(self.nodeCTs)
        if final_age == -1.:
            final_age = np.max(cts)
        node_ages = final_age * np.ones(cts.shape) - cts  # from creation time to age
        segs = self.get_segments()
        ages = np.zeros(segs.shape[0])
        for i, s in enumerate(segs):
            ages[i] = node_ages[s[1]]  # segment age based on second node (as in XylemFlux::linearSystem)
        return ages

    def get_nodes_index(self, ot):
        """ return node indices of segments with organ type @param ot """
        segments = self.get_segments()
        nodes = self.get_nodes()
        organTypes = self.get_organ_types()
        rootsegments = segments[organTypes == ot]
        rootsegments.flatten()
        np.sort(rootsegments, axis = None)
        nodesidx = np.unique(rootsegments)
        return nodesidx

    def get_nodes_organ_type(self, ot):
        """ return node coordinates of segments with organ type @param ot """
        nodes = self.get_nodes()
        return nodes[self.get_nodes_index(ot)]

    def get_segments_index(self, ot):
        """ return segment indices of organ type @param ot """
        organTypes = self.get_organ_types()
        segIdx = np.array(list(range(0, len(organTypes))))
        otsegs = segIdx[organTypes == ot]
        return otsegs
        
    def get_organ_types(self):
        """ segment organ types as numpy array """
        return np.array(self.organTypes)

    def get_organ_nodes_tips(self):
        """ return index of nodes at the end of each organ """
        organTypes = self.get_organ_types()
        segments = self.get_segments()
        get_y_node = lambda vec: vec[1]
        get_x_node = lambda vec: vec[0]
        get_nodetype = lambda y: organTypes[y - 1]
        nodesy = np.array([get_y_node(xi) for xi in segments], dtype = np.int64)
        nodesx = np.array([get_x_node(xi) for xi in segments], dtype = np.int64)
        nodesy = np.setdiff1d(nodesy, nodesx)  # select all the nodes which belong to tip of an organ
        nodes_type = np.array([get_nodetype(xi) for xi in nodesy], dtype = np.int64)
        tiproots = np.intersect1d(np.where(nodes_type == 2, nodesy, -1), nodesy)  # take root tips
        tipstem = np.intersect1d(np.where(nodes_type == 3, nodesy, -1), nodesy)  # take stem tips
        tipleaf = np.intersect1d(np.where(nodes_type == 4, nodesy, -1), nodesy)  # take leaf tips
        return tiproots, tipstem, tipleaf

    def get_organ_segments_tips(self):
        """ return index of segments at the end of each organ """
        tiproots, tipstems, tipleaves = self.get_organ_nodes_tips()
        tiproots = tiproots - np.ones(tiproots.shape, dtype = np.int64)  # segIndx = seg.y -1
        tipstems = tipstems - np.ones(tipstems.shape, dtype = np.int64)  # segIndx = seg.y -1
        tipleaves = tipleaves - np.ones(tipleaves.shape, dtype = np.int64)  # segIndx = seg.y -1
        return tiproots, tipstems, tipleaves

    def collar_index(self):
        """ returns the segment index of the collar segment """
        segs = self.segments
        for i, s in enumerate(segs):
            if s.x == 0:
                return i

    @staticmethod
    def read_rsml(file_name:str, verbose = True):
        """ reads an RSML file and converts to MappedSegments with units [cm]
        @file_name     the file name of the rsml, including file extension (e.g. "test.rsml" ) 
        @return a CPlantBox MappedSegments object
        """
        polylines, props, funcs, _ = rsml.read_rsml(file_name)
        bn = 0  # count base roots
        for i, _ in enumerate(polylines):
            if props["parent-poly"][i] < 0:
                bn += 1
        if bn > 1:
            rsml.artificial_shoot(polylines, props, funcs)
            if verbose:
                print("XylemFluxPython.read_rsml: added an artificial shoot")
        nodes, segs = rsml.get_segments(polylines, props)
        if verbose:
            print("XylemFluxPython.read_rsml: read rsml with", len(nodes), "nodes and", len(segs), "segments")
        nodes2 = [pb.Vector3d(n[0] , n[1] , n[2]) for n in nodes]  # Conversions to PlantBox types
        segs2 = [pb.Vector2i(int(s[0]), int(s[1])) for s in segs]

        radii, cts, types, tag_names = rsml.get_parameter(polylines, funcs, props)
        segRadii = np.zeros((segs.shape[0], 1))  # convert to paramter per segment
        segTypes = np.zeros((segs.shape[0], 1))
        for i, s in enumerate(segs):
            segRadii[i] = radii[s[1]]  # seg to node index
            segTypes[i] = types[s[1]]

        if verbose:
            print("                           cts [{:g}, {:g}] days".format(np.min(cts), np.max(cts)))
            print("                           raddii [{:g}, {:g}] cm".format(np.min(radii), np.max(radii)))
            print("                           subTypes [{:g}, {:g}] ".format(np.min(types), np.max(types)))
            print()

        return pb.MappedSegments(nodes2, cts, segs2, segRadii, segTypes)  # root system grid

    def get_incidence_matrix(self):
        """ returns the incidence matrix (number of segments)x(number of nodes) of the root system in self.ms 
        """
        segs = self.segments
        sn = len(segs)
        nn = len(self.nodes)  # TODO write getter
        ii_, jj_, vv_ = [], [], []
        for i, s in enumerate(segs):  # build incidence matrix from edges
            ii_.append(i)
            ii_.append(i)
            jj_.append(segs[i].x)
            jj_.append(segs[i].y)
            vv_.append(-1.)
            vv_.append(1.)
        return sparse.coo_matrix((np.array(vv_), (np.array(ii_), np.array(jj_))), shape = (sn, nn))

    def get_soil_matrix(self):
        """ maps nodes to soil matrix indices using the matrix B = (soil_matrix_indices) x (number_of_nodes-1); this is needed for upscaling  
               
            returns: 
            the matrix B
            a dictionary soil2matrix which maps soil_cell_index to soil_matrix_index
            a list matrix2soil which maps soil_matrix_index to soil_cell_index
        """
        soil2matrix = {}
        seg2cell = self.seg2cell
        segs = self.segments
        ns = len(segs)
        smi = 0  # soil matrix index
        ii_, jj_ = [], []

        for i, s in enumerate(segs):
            soil_cell_index = seg2cell[i]
            if not soil_cell_index in soil2matrix:
                soil2matrix[soil_cell_index] = smi
                smi += 1
            node_index = s.y
            ii_.append(soil2matrix[soil_cell_index])
            jj_.append(node_index - 1)

        B = sparse.coo_matrix((np.ones((len(ii_),)), (np.array(ii_), np.array(jj_))), shape = (smi, ns))

        cell_max = max(seg2cell.values()) + 1
        matrix2soil = np.zeros((smi,), dtype = np.int64)
        for i in range(0, cell_max):
            if i in soil2matrix:
                matrix2soil[soil2matrix[i]] = i

        return B, soil2matrix, matrix2soil


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

        seed = pb.Seed(self)

        # 2. Parse rsml, create self.static_organs
        
        parent_id = {}
        static_polys, static_types = [], []
        for i, root in enumerate(self.data.polylines):
            if types[i] in initial_sub_types:
                rp = self.getOrganRandomParameter(2, types[i])
                parent_id = self.data.properties["parent-poly"][i]
                pni = self.data.properties["parent-node"][i]

                length = self.polyline_length_(0, len(root) - 1, root)
                a = np.mean(radii[i])
                a_gr = rp.a_gr
                # print(a, len(radii[i]))
                r = 1.e6  # planted initially, and static
                lb = 0  # length of basal zone
                theta = 0  # insertion angle
                rlt = 1.e10  # root life time
                ln_ = []
                laterals = False
                if types[i] > 0:
                    p_survive = self.rand()
                    rlt_winter = rp.lambda_survive * ((-np.log(p_survive))**(1/rp.k_survive)) * 1225. # (25. + 10. * max(min(self.randn(),1.),-1.)) * 1225. #rp.lambda_survive * ((-np.log(rp.p_survive))**(1/rp.k_survive)) * 1225.
                else:
                    rlt_winter = 200. * 1225.
                print('rlt_winter',rlt_winter/1225)
                param = pb.RootSpecificParameter(types[i], lb, length, ln_, r, a, theta, rlt, laterals, a_gr, rlt_winter)  ############## which subType
                
                id = self.getOrganIndex()  # next index
                organ = pb.StaticRoot(id, param, length, pni)

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
                if parent_id < 0: # stem organ
                    assert types[i] == 0
                    seed.moveNode(n = pb.Vector3d(self.data.polylines[i][-1]), lId = 0)   
                    organ.addNode( n= pb.Vector3d(self.data.polylines[i][-1]),id = 0, t=0.)
                    for node in self.data.polylines[i][::-1][1:]: # reverse direction as stems go in the opposit direction
                        organ.addNode(n = pb.Vector3d(node),t =  0.)
                else:
                    # print(i, "parent", parent_id, "pni", pni)
                    parent = self.static_organs[parent_id]
                    organ.addNode(n = parent.getNode(pni), id = parent.getNodeId(pni), t = 0.)
                    # print(self.data.polylines[parent_id][pni])                  
                    for node in self.data.polylines[i]:
                        organ.addNode(n = pb.Vector3d(node), t = 0.)
                        
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
                    seed.addChild(organ)#self.static_organs[0])
        self.addOrgan(seed)

        # 4. The CPlantBox part
        # seed.initialize() # not called i.e. no tap root or basal roots are created
        self.oldNumberOfNodes = self.getNumberOfNodes()
        self.initCallbacks()
        
        self.nodes_py = np.array(list(map(lambda x: np.array(x), self.getNodes())))
        self.segments_py = np.array(list(map(lambda x: np.array(x), self.getSegments())), dtype = np.int64)

    def set_identical_laterals(self, initial_sub_types, lateral_subtypes, emerge_type):
        """ places laterals as in the original rsml, all start growing at once 
        
        initial_sub_types             subTypes of the initial static root system 
        lateral_subtypes              subTypes of the laterals wihtin the RSML
        emerge_type                   subType of the model lateral (todo: could be probabilistic)        
        """
        # where we add laterals on tope of the statics
        add_to_statics = np.array(initial_sub_types)[np.isin(initial_sub_types, lateral_subtypes)]
        ld = [] #lni, lt, , [], []
        ld1 = []
        types = self.data.properties[self.data.tagnames[2]]
        # print('types', types)
        for i, root in enumerate(self.data.polylines):
            if types[i] in lateral_subtypes:
                parent_id = self.data.properties["parent-poly"][i]
                # print('types[i]',types[i], 'types[parent_id]',types[parent_id], types[parent_id] in initial_sub_types)
                if types[parent_id] in initial_sub_types:
                    pni_init = self.data.properties["parent-node"][i] + 1 # why + 1?
                    parent = self.static_organs[parent_id]
                    pr = parent.getOrganRandomParameter()
                    init_num_kids = parent.getNumberOfChildren()  
                    #print('latId', parent.param().subType, pr.successorNo[0] , init_num_kids)
                    ##
                    # every time it sees a root in the rsml file, it will add pr.successorNo[0] laterals
                    ##
                    for latId in range(pr.successorNo[0] ):#(pr.successorNo[0] - init_num_kids)):
                        #print('in latid', latId)
                        pni = pni_init
                        creation_time = abs(self.randn() * pr.ldelays) # switch to rand?
                        # delay_mean = 0.
                        for kid_id in range(init_num_kids):
                            
                            kid = parent.getChild(kid_id)
                            if (kid.parentNI + 1 == pni) and (kid.getParameter('subType') == types[i]) and (np.isin(types[i], add_to_statics)):
                                # delay_mean = kid.getParameter('rlt_winter') 
                                creation_time = self.rand() * pr.ldelays
                                #creation_time = abs( (max(min(self.randn(),3.),-3.) / 3)* pr.ldelays) # delay_mean +
                                #print('delay1', creation_time, delay_mean,pr.ldelays,(max(min(self.randn(),3.),-3.) / 3),(max(min(self.randn(),3.),-3.) / 3))
                                ld1.append(creation_time)
                                pni -= 1
                                break
                                
                        #parent.getLatGrowthDelay()
                        # print('creation_time',creation_time,
                                # pr.successorOT,
                                # pr.successorST,
                                # pr.successorNo,
                                # pr.successorP,
                                # pr.successorP_age)
                        p_idx = 0#pr.getLateralType(pb.Vector3d(), 0, # ruleID: assumed to be 0 here
                        #                            creation_time)
                        #print('addlat','parent_st', parent.param().subType, 'kid_st' , pr.successorST[0][p_idx], 'pni',  pni)
                        if(p_idx >=0) :
                            emerge_type_ = pr.successorST[0][0][p_idx]
                            parent.addLateral(pni, emerge_type_, creation_time)
                            parent.param().laterals = True
                            # lt.append(emerge_type_)
                            # print("Root", i, ":", parent_id, pni, emerge_type, creation_time)
                        # else:
                        #    lt.append(p_idx)
                        ld.append(creation_time)
        
        self.nodes_py = np.array(list(map(lambda x: np.array(x), self.getNodes())))
        self.segments_py = np.array(list(map(lambda x: np.array(x), self.getSegments())), dtype = np.int64)                
        return ld, ld1

    def initialize_static_laterals(self):
        """ Creates lateral root instances from static roots lateralDelays, lateralTypes, and lateralDelays """
        for organ in self.static_organs.values():
            organ.initializeLaterals()
            # numKids = organ.getNumberOfChildren()
            # print('organ', organ.param().subType )
            # for numKid in range(numKids):
                # pr = organ.getChild(numKid).getOrganRandomParameter()
                # ps = organ.getChild(numKid).param()
                # print(ps.getK(),ps.subType)
            # raise Exception
        
        self.nodes_py = np.array(list(map(lambda x: np.array(x), self.getNodes())))
        self.segments_py = np.array(list(map(lambda x: np.array(x), self.getSegments())), dtype = np.int64)

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

