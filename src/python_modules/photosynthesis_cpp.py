import timeit


import numpy as np
#import matplotlib.pyplot as plt

import plantbox as pb
from plantbox import Photosynthesis
#import rsml_reader as rsml



class PhotosynthesisPython(Photosynthesis):
    """  wrapper for photosynthesis
       
    """

    def __init__(self, plant_, psiXylInit, ciInit):
        """ @param rs is either a pb.MappedRootSystem, pb.MappedSegments, or a string containing a rsml filename"""
        super().__init__( plant_, psiXylInit, ciInit)

    def get_nodes(self):
        """ converts the list of Vector3d to a 2D numpy array """
        return np.array(list(map(lambda x: np.array(x), self.rs.nodes)))

    def get_segments(self):
        """ converts the list of Vector2i to a 2D numpy array """
        return np.array(list(map(lambda x: np.array(x), self.rs.segments)), dtype = np.int64)

    def get_ages(self, final_age = 0.):
        """ converts the list of nodeCT to a numpy array of segment ages
        @param final_age [day]         current root system age, (default = 0 means detect from nodeCT)
        """
        cts = np.array(self.rs.nodeCTs)
        if final_age == 0.:
            final_age = np.max(cts)
        ages = final_age * np.ones(cts.shape) - cts  # from creation time to age
        return ages[1:]  # segment index is node index-1

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
        """ return node indices of organ type @param ot """
        organTypes = self.get_organ_types()
        segIdx = np.array(list(range(0, len(organTypes))))
        otsegs = segIdx[organTypes == ot]
        return otsegs

    def get_organ_types(self):
        """ segment organ types as numpy array """
        return np.array(self.rs.organTypes)

    def get_subtypes(self):
        """ segment sub types as numpy array """
        return np.array(self.rs.subTypes)

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

    def get_suf(self, sim_time):
        """ calculates the surface uptake fraction [1] of the root system at simulation time @param sim_time [day]
            (suf is constant for age independent conductivities)  """
        segs = self.rs.segments
        nodes = self.rs.nodes
        p_s = np.zeros((len(segs),))
        for i, s in enumerate(segs):
            p_s[i] = -500 - 0.5 * (nodes[s.x].z + nodes[s.y].z)  # constant total potential (hydraulic equilibrium)
        rx = self.solve_neumann(sim_time, -1.e5, p_s, False)  # False: matric potential not given per cell (but per segment), high number to recuce spurious fluxes
        print("rx ", np.min(rx), np.max(rx), np.mean(rx))
        fluxes = self.segFluxes(sim_time, rx, p_s, approx = False, cells = False)  # cm3/day, simTime,  rx,  sx,  approx, cells
        # print("fluxes ", np.min(fluxes) / -1.e5, np.max(fluxes) / -1.e5, np.mean(fluxes) / -1.e5)
        return np.array(fluxes) / -1.e5  # [1]

    def get_mean_suf_depth(self, sim_time):
        """  mean depth [cm] of water uptake based suf """
        suf = self.get_suf(sim_time)
        segs = self.rs.segments
        nodes = self.rs.nodes
        z_ = 0
        for i, s in enumerate(segs):
            z = 0.5 * (nodes[s.x].z + nodes[s.y].z)
            z_ += z * suf[i]
        return z_

    def find_base_segments(self):
        """ return all segment indices emerging from intial nodes given in self.dirichlet_ind
        (slow for large root systesms)
        """
        s_ = []
        segs = self.rs.segments
        for i, s in enumerate(segs):
            if s.x in self.dirichlet_ind:
                s_.append(i)
        print("XylemFluxPython.find_base_segments(): base segment indices for node indics", self.dirichlet_ind, "are", s_)
        return s_

    def get_krs(self, sim_time, seg_ind = [0]):
        """ calculatets root system conductivity [cm2/day] at simulation time @param sim_time [day] 
        if there is no single collar segment at index 0, pass indices using @param seg_ind, see find_base_segments        
        """
        segs = self.rs.segments
        nodes = self.rs.nodes
        p_s = np.zeros((len(segs),))
        for i, s in enumerate(segs):
            p_s[i] = -500 - 0.5 * (nodes[s.x].z + nodes[s.y].z)  # constant total potential (hydraulic equilibrium)
        rx = self.solve_dirichlet(sim_time, -15000, 0., p_s, cells = False)
        jc = 0
        for i in seg_ind:
            jc -= self.axial_flux(i, sim_time, rx, p_s, [], cells = False, ij = True)
        krs = jc / (-500 - 0.5 * (nodes[segs[0].x].z + nodes[segs[0].y].z) - rx[0])
        return krs , jc

    def get_eswp(self, sim_time, p_s):
        """ calculates the equivalent soil water potential [cm] at simulation time @param sim_time [day] for 
        the soil matric potential @param p_s [cm] given per cell """
        segs = self.rs.segments
        nodes = self.rs.nodes
        seg2cell = self.rs.seg2cell
        suf = self.get_suf(sim_time)
        eswp = 0.
        for i, s in enumerate(segs):
            eswp += suf[i] * (p_s[seg2cell[i]] + 0.5 * (nodes[s.x].z + nodes[s.y].z))  # matric potential to total potential
        return eswp

    def kr_f(self, age, st, ot = 2 , numleaf = 2, seg_ind = 0):
        """ root radial conductivity [1 day-1] for backwards compatibility """
        return self.kr_f_cpp(seg_ind, age, st, ot, numleaf)  # kr_f_cpp is XylemFlux::kr_f

    def kx_f(self, age, st, ot = 2, seg_ind = 0):
        """ root axial conductivity [cm3 day-1]  for backwards compatibility """
        return self.kx_f_cpp(seg_ind, age, st, ot)  # kx_f_cpp is XylemFlux::kx_f

    @staticmethod
    def convert_(x, dtype = np.float64):
        """ not used anymore (?) """
        return np.array(list(map(lambda x: np.array(x, dtype), x)), dtype)  # is there a better way?

