import timeit


import numpy as np
#import matplotlib.pyplot as plt

import plantbox as pb
from plantbox import PhloemFlux
from functional.xylem_flux import XylemFluxPython
#import rsml_reader as rsml



class PhloemFluxPython(PhloemFlux):
    """  wrapper for photosynthesis
       somehow, double inheritance is not working? so i m copying the 
       necessary XylemFluxPython functions
    """

    def __init__(self, plant_, psiXylInit, ciInit):
        """ @param rs is either a pb.MappedRootSystem, pb.MappedSegments, or a string containing a rsml filename"""
        super().__init__( plant_, psiXylInit, ciInit)
    
    def getPsiAir(self,RH, TairC):#constants are within photosynthesys.h
        return np.log(RH) * self.rho_h2o * self.R_ph * (TairC + 237.3)/self.Mh2o * (1/0.9806806)  ; #in cm
        
    def axial_flux(self, seg_ind, sim_time, rx, sxx, k_soil = [], cells = True, ij = True):
        """ returns the exact axial flux of segment ij of xylem model solution @param rx
            @param seg_ind              segment index 
            @param sim_time [day]       needed for age dependent conductivities (age = sim_time - segment creation time)        
            @param rx [cm]              root xylem matric potentials per root system node
            @param sxx [cm]             soil matric potentials given per segment or per soil cell
            @param k_soil [day-1]       optionally, soil conductivities can be prescribed per segment, 
                                        conductivity at the root surface will be limited by the value, i.e. kr = min(kr_root, k_soil) 
            @param cells                indicates if the matric potentials are given per cell (True) or by segments (False)
            @param ij                   True: calculate axial flux in node i, False: in node j; note that they are not equal due to radial fluxes 
            @return [cm3 day-1] axial volumetric flow rate             
        """
        s = self.rs.segments[seg_ind]
        nodes = self.rs.nodes
        numleaf = 0
        organTypes = self.get_organ_types()
        ot = int(organTypes[seg_ind])  # for conductivities kr, kx
        assert ot == 2 #this function only works for root (and stem?) segments
        if len(k_soil) > 0:
            ksoil = k_soil[seg_ind]
        else:
            ksoil = 1.e9  # much
        # node x and y of stem and leave segments are revers with regards to the nodes x and y of roots.
        if sum(((not ij) , (ot == 4) or (ot == 3))) == 1:  # check if only one of the condition is true
            j, i = int(s.x), int(s.y)  # node indices
        else:
            i, j = int(s.x), int(s.y)
        n1, n2 = self.rs.nodes[i], self.rs.nodes[j]  # nodes
        v = n2.minus(n1)
        l = v.length()  # length of segment
        v.normalize()  # normalized v.z is needed for qz
        if cells:
            cell_ind = self.rs.seg2cell[seg_ind]
            if cell_ind >= 0:  # y node belowground
                if len(sxx) > 1:
                    p_s = sxx[cell_ind]  # soil pressure at collar segment
                else:
                    p_s = sxx[0]
            else:
                p_s = self.airPressure
        else:
            p_s = sxx[seg_ind]

        a = self.rs.radii[seg_ind]  # radius
        st = int(self.rs.subTypes[seg_ind])  # conductivities kr, kx
        age = sim_time - self.rs.nodeCTs[int(s.y)]
        if ot == 4:  # to know which x-th leaf segment it is, to fetch the right gs value
                indices = [i for i, x in enumerate(organTypes) if x == 4]
                numleaf = indices.index(seg_ind)
                if self.pg[0] != 0:
                    p_s = self.pg[numleaf]
                    
        kr = self.kr_f(age, st, ot, seg_ind, cells)  # c++ conductivity call back functions
        kr = min(kr, ksoil)
        kx = self.kx_f(age, st, ot,seg_ind)
        
        if a * kr > 1.e-16:
            tau = np.sqrt(2 * a * np.pi * kr / kx)  # cm-2
            AA = np.array([[1, 1], [np.exp(tau * l), np.exp(-tau * l)] ])
            bb = np.array([rx[i] - p_s, rx[j] - p_s])  # solve for solution
            d = np.linalg.solve(AA, bb)  # compute constants d_1 and d_2 from bc
            dpdz0 = d[0] * tau - d[1] * tau  # insert z = 0, z = l into exact solution
        else:  # solution for kr = 0, or a = 0
            dpdz0 = (rx[j] - rx[i]) / l

        f = kx * (dpdz0 + v.z)
        if ij:
            f = f * (-1)
        print(f)
        return f
        
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
        
    def test(self):
        """ perfoms some sanity checks, and prints to the console """
        print("\nPhloemFluxPython.test():")
        # 1 check if segment index is node index-1
        segments = self.get_segments()
        nodes = self.get_nodes()
        for i, s_ in enumerate(segments):
            if i != s_[1] - 1:
                raise "Segment indices are mixed up"
        print(len(segments), "segments")
        # 2 check for very small segments
        seg_length = self.rs.segLength()
        c = 0
        for i, l in enumerate(seg_length):
            if l < 1e-5:
                print(i, l)
                c += 1
        print(c, "segments with length < 1.e-5 cm")
        # 3 check for type range, index should start at 0
        types = self.rs.subTypes
        if np.min(types) > 0:
            print("Warning: types start with index", np.min(types), "> 0 !")
        print("{:g} different root types from {:g} to {:g}".format(np.max(types) - np.min(types) + 1, np.min(types), np.max(types)))
        # 4 Print segment age range
        ages = self.get_ages()
        print("ages from {:g} to {:g}".format(np.min(ages), np.max(ages)))
        # 4 check for unmapped indices
        map = self.rs.seg2cell
        organTypes = self.rs.organTypes
        for seg_id, cell_id in map.items():
            if (cell_id < 0) and organTypes[seg_id] == 2:
                print("Warning: root segment ", seg_id,"is NOT mapped, this may cause problems with coupling!", nodes[segments[seg_id][0]], nodes[segments[seg_id][1]])
            if (cell_id >= 0) and organTypes[seg_id] != 2:
                print("Warning: shoot segment ", seg_id, "organ type",organTypes[seg_id],"IS mapped, this may cause problems with coupling!", nodes[segments[seg_id][0]], nodes[segments[seg_id][1]])
                

    def kr_f(self, age, st, ot = 2 , seg_ind = 0, cells = False):
        """ root radial conductivity [1 day-1] for backwards compatibility """
        return self.kr_f_cpp(seg_ind, age, st, ot, cells)  # kr_f_cpp is XylemFlux::kr_f

    def kx_f(self, age, st, ot = 2, seg_ind = 0):
        """ root axial conductivity [cm3 day-1]  for backwards compatibility """
        return self.kx_f_cpp(seg_ind, age, st, ot)  # kx_f_cpp is XylemFlux::kx_f

    @staticmethod
    def convert_(x, dtype = np.float64):
        """ not used anymore (?) """
        return np.array(list(map(lambda x: np.array(x, dtype), x)), dtype)  # is there a better way?

