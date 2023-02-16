import timeit

import numpy as np
from scipy import sparse
import scipy.sparse.linalg as LA
import matplotlib.pyplot as plt

import plantbox as pb
from plantbox import XylemFlux
from xylem_flux import XylemFluxPython


class HydraulicsDoussan(XylemFluxPython):
    """  Root hydraulic model (following Doussn et al. )
    
    TODO
    
    """

    def __init__(self, rs):
        """ @param rs is either a pb.MappedRootSystem, pb.MappedSegments, or a string containing a rsml filename"""
        super().__init__(rs)

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

    def doussan_system(self, sim_time, sxx, cells = True, soil_k = []):
        """ sets up the system matrix self.A and load vector"self.b """
        IM = self.get_incidence_matrix()
        IMt = IM.transpose()
        kx_ = np.divide(self.getKx(sim_time), self.rs.segLength())  # / dl
        Kx = sparse.diags(kx_)
        kr = np.zeros((IM.shape[1],))
        kr[1:] = np.array(self.getEffKr(sim_time))  #  change to getKr(), use that instead
        Kr = sparse.diags(kr)
        A = IMt @ Kx @ IM + Kr  # Laplacian IMt @ Kx @ IM =: L in Hess paper
        return A

        # if cells:
        #     hs = np.zeros((IM.shape[1],))
        #     hs[1:] = self.getHs(sxx)
        #     self.b = Kr * hs
        # else:
        #     hs = np.zeros((IM.shape[1],))
        #     hs[1:] = sxx
        #     self.b = Kr * hs

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

    def axial_flux(self, seg_ind, sim_time, rx):
        """ returns the exact axial flux of segment ij of xylem model solution @param rx
            @param seg_ind              segment index 
            @param sim_time [day]       needed for age dependent conductivities (age = sim_time - segment creation time)        
            @param rx [cm]              root xylem matric potentials per root system node
            @return [cm3 day-1] axial volumetric flow rate             
        """
        s = self.rs.segments[seg_ind]
        i, j = s.x, s.y
        a = self.rs.radii[seg_ind]  # radius
        st = int(self.rs.subTypes[seg_ind])  # sub type
        age = sim_time - self.rs.nodeCTs[int(s.y)]
        kr = self.kr_f(age, st)  # c++ conductivity call back functions
        kx = self.kx_f(age, st)  # c++ conductivity call back functi
        dpdz0 = (rx[j] - rx[i]) / l
        f = -kx * (dpdz0)
        return f

