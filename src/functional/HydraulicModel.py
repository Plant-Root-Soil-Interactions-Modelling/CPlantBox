import timeit

import numpy as np
from scipy import sparse
import scipy.sparse.linalg as LA
import matplotlib.pyplot as plt

import plantbox as pb
from functional.xylem_flux import XylemFluxPython as XylemFluxPythonCpp

import rsml.rsml_reader as rsml


class PlantHydraulicModel():  # PlantHydraulicModelCpp
    """  Root hydraulic models classs 
    
    """

    def __init__(self, method, ms):
        """ @param rs is  pb.MappedSegments, or a string containing a rsml filename"""
        if isinstance(rs, str):
            rs = self.read_rsml(rs)
            super().__init__(rs)
        else:
            super().__init__(rs)

        if method == "Meunier":
            pass
        elif method == "Doussan":
            pass
        else:
            raise ValueError("PlantHydraulicModel(): Unknown method: " + method)

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
        return pb.MappedSegments(nodes2, cts, segs2, segRadii, segTypes)  # root system grid

    def get_incidence_matrix(self):
        """ retruns the incidence matrix (number of segments, number of nodes) of the root system in self.rs """
        segs = self.rs.segments
        sn = len(segs)
        nn = len(self.rs.nodes)  # TODO write getter
        ii_, jj_, vv_ = [], [], []
        for i, s in enumerate(segs):  # build incidence matrix from edges
            ii_.append(i)
            ii_.append(i)
            jj_.append(segs[i].x)
            jj_.append(segs[i].y)
            vv_.append(-1.)
            vv_.append(1.)
        return sparse.coo_matrix((np.array(vv_), (np.array(ii_), np.array(jj_))), shape = (sn, nn))

    def doussan_system_matrix(self, sim_time):
        """ 
            returns all that is needed to solve root hydraulics with Doussan method:             
            sparse System matrix A,
            Kr diag kr values
            K_x,collar        
            A.dot(H_x) = Kr.dot(H_{s,r}) + K_x,collar * e_collar * wilting_point
        """
        IM = self.get_incidence_matrix()
        IMt = IM.transpose()
        kx_ = np.divide(self.getKx(sim_time), self.rs.segLength())  # / dl
        Kx = sparse.diags(kx_)
        kr = np.array(self.getEffKr(sim_time))
        kr = np.maximum(np.ones(kr.shape) * 1.e-12, kr)  #  limit to a small value for inversion
        Kr = sparse.diags(kr)
        L = IMt @ Kx @ IM  # Laplacian
        L_ = L[1:, 1:].tocsc()
        return  L_ + Kr, Kr, kx_[self.collar_index()]  # L_{N-1} + Kr, se Hess paper

    def get_soil_matrix(self):
        """ maps nodes to soil matrix indices using the matrix B = (soil_matrix_indices) x (number_of_nodes-1) 
               
            returns: 
            the matrix B
            a dictionary soil2matrix which maps soil_cell_index to soil_matrix_index
            a list matrix2soil which maps soil_matrix_index to soil_cell_index
        """
        soil2matrix = {}
        seg2cell = self.rs.seg2cell
        segs = self.rs.segments
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

    def collar_index(self):
        """ returns the segment index of the collar segment, node index of the collar node is 0 """
        segs = self.rs.segments
        for i, s in enumerate(segs):
            if s.x == 0:
                return i








    def get_transpiration(self):
        
        
    # transpiration  # sum of radial fluxes
    # radial_fluxes
    # axial_fluxes
    # xylem_potentials

    def axial_flux(self, seg_ind, sim_time, rx):
        """ returns the exact axial flux of segment ij of xylem model solution @param rx
            @param seg_ind              segment index 
            @param sim_time [day]       needed for age dependent conductivities (age = sim_time - segment creation time)        
            @param rx [cm]              root xylem matric potentials per root system node
            @return [cm3 day-1] axial volumetric flow rate             
        """
        s = self.rs.segments[seg_ind]
        i, j = s.x, s.y
        n1, n2 = self.rs.nodes[i], self.rs.nodes[j]  # nodes
        v = n2.minus(n1)
        l = v.length()
        a = self.rs.radii[seg_ind]  # radius
        st = int(self.rs.subTypes[seg_ind])  # sub type
        age = sim_time - self.rs.nodeCTs[int(s.y)]
        kr = self.kr_f(age, st)  # c++ conductivity call back functions
        kx = self.kx_f(age, st)  # c++ conductivity call back functi
        dpdz0 = (rx[j] - rx[i]) / l
        print(seg_ind, rx[j], rx[i], s, l, kx)
        f = -kx * (dpdz0)
        return f

    def collar_flux(self, sim_time, rx):
        """ returns the exact transpirational flux of the xylem model solution @param rx
            @see axial_flux        
        """
        return self.axial_flux(self.collar_index(), sim_time, rx)

    def axial_fluxes(self, sim_time, rx):
        """ returns the axial fluxes 
        @see axial_flux  
        """
        n = len(self.rs.segments)  # TODO getter
        return np.array([self.axial_flux(i, sim_time, rx) for i in range(0, n)])

    def get_krs(self, sim_time):
        """ calculatets root system conductivity [cm2/day] at simulation time @param sim_time [day]             
        """
        n = len(self.rs.segments)  # TODO getter
        H_sr = np.ones((n,)) * (-500)

        A, Kr, Kx0 = self.doussan_system_matrix(sim_time)  # A.dot(H_x) = Kr.dot(H_{s,r}) + K_x,collar * e_collar * wilting_point
        b = Kr.dot(H_sr)
        # print("b", b.shape)
        b[self.collar_index()] += Kx0 * (-15000)
        # print("collar_index", self.collar_index())
        Hx = LA.spsolve(A, b)

        # print("Hx", np.min(Hx), np.max(Hx))

        seg_ind = self.collar_index()
        s = self.rs.segments[seg_ind]
        i, j = s.x, s.y
        n1, n2 = self.rs.nodes[i], self.rs.nodes[j]  # nodes
        v = n2.minus(n1)
        l = v.length()
        a = self.rs.radii[seg_ind]  # radius
        st = int(self.rs.subTypes[seg_ind])  # sub type
        age = sim_time - self.rs.nodeCTs[int(s.y)]
        kx = self.kx_f(age, st)
        dpdz0 = (Hx[seg_ind] - (-15000)) / l
        print(seg_ind, Hx[seg_ind], s, a, l, kx)

        t_act = -kx * (dpdz0)
        # t_act = self.collar_flux(sim_time, Hx)

        # print("t_act", t_act, "Hx0", Hx[0])
        krs = -t_act / ((-500) - Hx[seg_ind])
        return krs, t_act

    def get_krs_old(self, sim_time, seg_ind = [0]):
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
            jc -= self.axial_flux(i, sim_time, rx)

        krs = jc / (-500 - 0.5 * (nodes[segs[0].x].z + nodes[segs[0].y].z) - rx[self.dirichlet_ind[0]])
        return krs , jc

