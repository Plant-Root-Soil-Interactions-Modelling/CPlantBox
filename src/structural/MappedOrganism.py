"""
Wrapper for MappedPlant, MappedSegments, or MappedOrganism 
to add functionality that is easier implemented in Python 
"""
import numpy as np
from scipy import sparse


class MappedPlantPython():

    def __init__(self, ms):
        """ @param ms (MappedPlant, MappedSegments, or MappedOrganism) to be wrapped """
        self.ms = ms

    def toNumpy(self,cpbArray):
        """ converts the cpbArray to a numpy array """
        return np.array(list(map(lambda x: np.array(x), cpbArray)))
        
    def get_nodes(self):
        """ converts the list of Vector3d to a 2D numpy array """
        return np.array(list(map(lambda x: np.array(x), self.ms.nodes)))

    def get_segments(self):
        """ converts the list of Vector2i to a 2D numpy array """
        return np.array(list(map(lambda x: np.array(x), self.ms.segments)), dtype = np.int64)

    def get_subtypes(self):
        """ segment sub types as numpy array """
        return np.array(self.ms.subTypes)

    def get_ages(self, final_age = -1.):
        """ converts the list of nodeCT to a numpy array of segment ages
        @param final_age [day]         current root system age, (default = 0 means detect maximum from nodeCT)
        """
        cts = np.array(self.ms.nodeCTs)
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
        return np.array(self.ms.organTypes)

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
        segs = self.ms.segments
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
        segs = self.ms.segments
        sn = len(segs)
        nn = self.ms.getNumberOfMappedNodes()  
        # ii_, jj_, vv_ = [], [], []
        # for i, s in enumerate(segs):  # build incidence matrix from edges
        #     ii_.append(i)
        #     ii_.append(i)
        #     jj_.append(segs[i].x)
        #     jj_.append(segs[i].y)
        #     vv_.append(-1.)
        #     vv_.append(1.)
        # # first optimization 
        # ii_ = np.repeat(np.arange(sn), 2)
        # jj_ = np.array([(s.x, s.y) for s in segs]).ravel()
        # vv_ = np.tile([-1., 1.], sn) 
        # # second optimization
        ii_ = np.repeat(np.arange(sn, dtype=np.int32), 2)
        # Build jj_ efficiently using fromiter
        xs = np.fromiter((s.x for s in segs), dtype=np.int32, count=sn) # TODO C++ gettter could improve performance
        ys = np.fromiter((s.y for s in segs), dtype=np.int32, count=sn) 
        jj_ = np.empty(sn * 2, dtype=np.int32)
        jj_[0::2] = xs
        jj_[1::2] = ys
        # Values array
        vv_ = np.empty(sn * 2, dtype=np.float64)
        vv_[0::2] = -1.0
        vv_[1::2] =  1.0        
        return sparse.coo_matrix((np.array(vv_), (np.array(ii_), np.array(jj_))), shape = (sn, nn))

    def get_soil_matrix(self):
        """ maps nodes to soil matrix indices using the matrix B = (soil_matrix_indices) x (number_of_nodes-1); this is needed for upscaling  
               
            returns: 
            the matrix B
            a dictionary soil2matrix which maps soil_cell_index to soil_matrix_index
            a list matrix2soil which maps soil_matrix_index to soil_cell_index
        """
        soil2matrix = {}
        seg2cell = self.ms.seg2cell
        segs = self.ms.segments
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

