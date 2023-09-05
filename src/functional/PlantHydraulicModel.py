import timeit

import numpy as np
from scipy import sparse
import scipy.sparse.linalg as LA
import matplotlib.pyplot as plt

import plantbox as pb
from plantbox import PlantHydraulicModel as PlantHydraulicModelCPP

import rsml.rsml_reader as rsml


class PlantHydraulicModel(PlantHydraulicModelCPP):
    """  Root hydraulic models classs 
    
    
        will incorporate Doussan and Meunier (currently working on Doussan)
    """

    def __init__(self, method, rs, params):
        """ 
        @param method is either "Meunier" or "Doussan" 
        @param ms is of type MappedSegments, or a string containing a rsml filename
        @param params hydraulic conductivities described by PlantHydraulicParameters
        """
        if isinstance(rs, str):
            rs = self.read_rsml(rs)
            super().__init__(rs, params)
        else:
            super().__init__(rs, params)

        self.method = method
        if not method in ["Doussan"]:  # "Meunier",
            raise ValueError("PlantHydraulicModel(): Unknown method: " + method)

    @staticmethod
    def read_rsml(file_name:str, verbose = True):
        """ 
        reads an RSML file and converts to MappedSegments with units [cm]
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

    def get_nodes(self):
        """ converts the list of Vector3d to a 2D numpy array """
        return np.array(list(map(lambda x: np.array(x), self.rs.nodes)))

    def get_segments(self):
        """ converts the list of Vector2i to a 2D numpy array """
        return np.array(list(map(lambda x: np.array(x), self.rs.segments)), dtype = np.int64)

    def get_efffective_kr(self, sim_time):
        """ effective radial conductivities per segment [cm2 day-1] """
        return np.array(self.params.getEff(self.rs, sim_time))

    def get_kr(self, sim_time):
        """ radial conductivities per segment [1 day-1] """
        return np.array(self.params.getKr(self.rs, sim_time))

    def get_kx(self, sim_time):
        """ axial conductivities per segment [cm3 day-1]"""
        return np.array(self.params.getKx(self.rs, sim_time))

    def get_hs(self, sx):
        """ soil matric potential per segment [cm] """
        return np.array(self.params.getHs(self.rs, sx))

    def get_incidence_matrix(self):
        """ returns the incidence matrix (number of segments, number of nodes) of the root system in self.rs 
            
            TODO move to MappedSegments 
        """
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

    def get_soil_matrix(self):
        """ maps nodes to soil matrix indices using the matrix B = (soil_matrix_indices) x (number_of_nodes-1), (needed for upscaling)
               
            returns: 
            the matrix B
            a dictionary soil2matrix which maps soil_cell_index to soil_matrix_index
            a list matrix2soil which maps soil_matrix_index to soil_cell_index
            
            TODO move to mapped segments?
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

    def get_doussan_system(self, sim_time):
        """ 
            returns all that is needed to solve root hydraulics with Doussan method TODO ref paper...             
            sparse System matrix A,
            Kr diag kr values
            K_x,collar        
            A.dot(H_x) = Kr.dot(H_{s,r}) + K_x,collar * e_collar * wilting_point
        """
        IM = self.get_incidence_matrix()
        IMt = IM.transpose()
        kx_ = np.divide(self.params.getKx(self.rs, sim_time), self.rs.segLength())  # / dl
        Kx = sparse.diags(kx_)
        kr = np.array(self.params.getEffKr(self.rs, sim_time))
        kr = np.maximum(np.ones(kr.shape) * 1.e-12, kr)  #  limit to a small value for inversion
        Kr = sparse.diags(kr)
        L = IMt @ Kx @ IM  # Laplacian
        L_ = L[1:, 1:].tocsc()
        return  L_ + Kr, Kr, kx_[self.collar_index()]  # L_{N-1} + Kr

    def collar_index(self):
        """ returns the segment index of the collar segment, node index of the collar node is 0 """
        segs = self.rs.segments
        for i, s in enumerate(segs):
            if s.x == 0:
                return i

    def update(self, sim_time):  # rs_age + simtime...
        """ call before solve() """
        A_d, self.Kr, self.kx0 = self.get_doussan_system(sim_time)
        self.ci = self.collar_index()
        A_n = A_d.copy()
        A_n[self.ci, self.ci] -= self.kx0
        print("update(): invert matrix start (splu)")
        self.A_n_splu = LA.splu(A_n)
        self.A_d_splu = LA.splu(A_d)
        print("update(): invert done.")

    def solve(self, rsx, t_pot, wilting_point):
        """ solves the hydraulic model
            @param t_pot [cm3 day-1]     potential transpiration rate
            @param rsx [cm]              soil total potential around root collar, if it is below the wilting_point, 
                                         dirichlet boundary conditions are assumed. Set sx = 0 to disable this behaviour.
            @parm wiltingPoint [cm]      the plant wilting point   
            @return [cm] root xylem total potential
        """
        b = self.Kr.dot(rsx)
        b[self.ci, 0] += self.kx0 * wilting_point
        # rx = np.expand_dims(sparse.linalg.spsolve(A_d, b), axis = 1)
        rx = self.A_d_splu.solve(b)
        q_dirichlet = -self.Kr.dot(rsx - rx)  # both total potentials
        if np.sum(q_dirichlet) <= t_pot:
            b = self.Kr.dot(rsx)
            b[self.ci, 0] += t_pot
            # rx = np.expand_dims(sparse.linalg.spsolve(A_n, b), axis = 1)
            rx = self.A_n_splu.solve(b)
        return rx

    def radial_fluxes(self, rx, rsx):
        """ returns the radial fluxes [cm3 day-1]"""
        return -self.Kr.dot(rsx - rx)

    def get_transpiration(self, rx, rsx):
        """ actual transpiration [cm3 day-1]"""
        return np.sum(self.radial_fluxes(rx, rsx))

    def get_krs(self, sim_time):
        """ calculatets root system conductivity [cm2/day] at simulation time @param sim_time [day]             
        """
        n = len(self.rs.segments)  # TODO getter
        s = self.rs.segments[self.ci]
        n2 = self.rs.nodes[s.y]
        rsx = np.ones((n, 1)) * (-500)
        b = self.Kr.dot(rsx)
        b[self.ci, 0] += self.kx0 * -15000
        rx = self.A_d_splu.solve(b)  # total matric potential
        t_act = self.get_transpiration(rx, rsx)
        krs = -t_act / ((-500) - (rx[self.ci, 0] - n2.z))  # from total to matric
        return krs, t_act

    def get_suf(self):
        """ standard uptake fraction per root segment, 
        call update(sim_time) before using it """
        n = len(self.rs.segments)  # TODO getter
        rsx = np.ones((n, 1)) * (-500)
        b = self.Kr.dot(rsx)
        b[self.ci, 0] += self.kx0 * -15000
        rx = self.A_d_splu.solve(b)
        q = self.radial_fluxes(rx, rsx)
        return np.array(q) / np.sum(q)
