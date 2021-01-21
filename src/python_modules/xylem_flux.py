import timeit
import math

import numpy as np
from scipy import sparse
import scipy.sparse.linalg as LA

import plantbox as pb
from plantbox import XylemFlux
import rsml_reader as rsml  # todo


class XylemFluxPython(XylemFlux):
    """  Hybrid flux solver (following Meunier et al.)
    
        C++ part is defined in CPlantBox XylemFlux.hh and XylemFlux.cpp
        
        Calculates water movement within the xylems, assuming a constant matric potential around each xylem segment,
        given for each segment, or for each soil cell. 
        
        The root surface flux is calculated exactly as in Meunier et al.        
    """

    def __init__(self, rs):
        """ @param rs is either a pb.MappedRootSystem or a string containing a rsml filename"""
        if isinstance(rs, str):
            rs = self.read_rsml(rs)
            super().__init__(rs)
        else:
            super().__init__(rs)
        
        self.seg_ind = [0]  # segment indices for Neuman flux
        self.node_ind = [0]  # node indices for Dirichlet flux

    def solve_neumann(self, sim_time :float, value, sxx, cells :bool, soil_k=[]) :
        """ solves the flux equations, with a neumann boundary condtion, see solve()
            @param sim_time [day]       needed for age dependent conductivities (age = sim_time - segment creation time)
            @param value [cm3 day-1]    tranpirational flux is negative
            @param sxx [cm]             soil matric potentials given per segment or per soil cell
            @param cells                indicates if the matric potentials are given per cell (True) or by segments (False)
            @param soil_k [day-1]       optionally, soil conductivities can be prescribed per segment, 
                                        conductivity at the root surface will be limited by the value, i.e. kr = min(kr_root, k_soil)  
            @return [cm] root xylem pressure per root system node         
         """
        # start = timeit.default_timer()
        if isinstance(value, (float, int)):
            n = len(self.seg_ind)
            value = [value / n] * n
            
        if len(soil_k) > 0:
            self.linearSystem(sim_time, sxx, cells, soil_k)  # C++ (see XylemFlux.cpp)
        else:
            self.linearSystem(sim_time, sxx, cells)  # C++ (see XylemFlux.cpp)
        Q = sparse.coo_matrix((np.array(self.aV), (np.array(self.aI), np.array(self.aJ))))
        Q = sparse.csr_matrix(Q)
        Q, b = self.bc_neumann(Q, self.aB, self.seg_ind, value)  # cm3 day-1
        x = LA.spsolve(Q, b, use_umfpack=True)  # direct
        # print ("linear system assembled and solved in", timeit.default_timer() - start, " s")
        return x

    def solve_dirichlet(self, sim_time :float, value :list, sxc :float, sxx, cells :bool, soil_k=[]):
        """ solves the flux equations, with a dirichlet boundary condtion, see solve()
            @param sim_time [day]     needed for age dependent conductivities (age = sim_time - segment creation time)
            @param scx                depricated (unused)
            @param value [cm]         root collar pressure head 
            @param sxx [cm]           soil matric potentials given per segment or per soil cell
            @param cells              indicates if the matric potentials are given per cell (True) or by segments (False)
            @param soil_k [day-1]     optionally, soil conductivities can be prescribed per segment, 
                                      conductivity at the root surface will be limited by the value, i.e. kr = min(kr_root, k_soil)  
            @return [cm] root xylem pressure per root system node
         """
        if isinstance(value, (float, int)):
            n = len(self.node_ind)
            value = [value] * n       
         
        if len(soil_k) > 0:
            self.linearSystem(sim_time, sxx, cells, soil_k)  # C++ (see XylemFlux.cpp)
        else:
            self.linearSystem(sim_time, sxx, cells)
        self.linearSystem(sim_time, sxx, cells)  # C++ (see XylemFlux.cpp)
        Q = sparse.coo_matrix((np.array(self.aV), (np.array(self.aI), np.array(self.aJ))))
        Q = sparse.csr_matrix(Q)
        Q, b = self.bc_dirichlet(Q, self.aB, self.node_ind, value)
        x = LA.spsolve(Q, b, use_umfpack=True)
        return x

    def solve(self, sim_time :float, trans :list, sx :float, sxx, cells :bool, wilting_point :float, soil_k=[]):
        """ solves the flux equations using Neumann and switching to dirichlet in case wilting point is reached in root collar 
            @param sim_time [day]        needed for age dependent conductivities (age = sim_time - segment creation time)
            @param trans [cm3 day-1]     transpiration rate
            @param sx [cm]               soil matric potential around root collar, if it is below the wilting_point, 
                                         dirichlet boundary conditions are assumed. Set sx = 0 to disable this behaviour.
            @param sxx [cm]              soil matric potentials given per segment or per soil cell
            @param cells                 indicates if the matric potentials are given per cell (True) or by segments (False)
            @parm wiltingPoint [cm]      the plant wilting point   
            @param soil_k [day-1]       optionally, soil conductivities can be prescribed per segment, 
                                        conductivity at the root surface will be limited by the value, i.e. kr = min(kr_root, k_soil)  
            @return [cm] root xylem pressure per root system node
        """
        eps = 1

        if sx >= wilting_point - eps:

            x = self.solve_neumann(sim_time, trans, sxx, cells, soil_k)  # try neumann, if below wilting point, switch to Dirichlet

            if x[0] <= wilting_point:
                Q = sparse.coo_matrix((np.array(self.aV), (np.array(self.aI), np.array(self.aJ))))
                Q = sparse.csr_matrix(Q)
                Q, b = self.bc_dirichlet(Q, self.aB, [0], [float(wilting_point)])
                try:
                    x = LA.spsolve(Q, b, use_umfpack=True)
                except:
                    print("XylemFluxPython.solve: Exeption solving Dirichlet")
                    print("Dirichlet at ", trans - sx, "cm")
                    print("b", b)
                    raise
        else:
            print("XylemFluxPython.solve: used Dirichlet because collar cell soil matric potential is below wilting point", sx)
            x = self.solve_dirichlet(sim_time, wilting_point, sx, sxx, cells, soil_k)

        return x

    def axial_flux(self, seg_ind, sim_time, rx, sxx, k_soil=[], cells=True, ij=True):
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
        organTypes = self.get_organ_types()  # collar segment
        ot = int(organTypes[seg_ind])  # conductivities kr, kx
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
        if ot == 4 and len(self.gs) > 0:  # to know which x-th leaf segment it is, to fetch the right gs value 
                indices = [i for i, x in enumerate(organTypes) if x == 4]
                numleaf = indices.index(seg_ind)
        kr = self.kr_f(age, st, ot, numleaf)  # c++ conductivity call back functions
        kr = min(kr, ksoil)
        kx = self.kx_f(age, st, ot)
        tau = math.sqrt(2 * a * math.pi * kr / kx)  # cm-2
        AA = np.array([[1, 1], [math.exp(tau * l), math.exp(-tau * l)] ])          
        bb = np.array([rx[i] - p_s, rx[j] - p_s])  # solve for solution
        d = np.linalg.solve(AA, bb)  # compute constants d_1 and d_2 from bc
        dpdz0 = d[0] * tau - d[1] * tau  # insert z = 0, z = l into exact solution   
        f = kx * (dpdz0 + v.z)
        if ij :
            f = f * (-1)
        return f  

    def collar_flux(self, sim_time, rx, sxx, k_soil=[], cells=True):
        """ returns the exact transpirational flux of the xylem model solution @param rx
            @see axial_flux        
        """
        return self.axial_flux(0, sim_time, rx, sxx, k_soil, cells, True)

    def axial_fluxes(self, sim_time, rx, sxx, k_soil=[], cells=True):
        """ returns the axial fluxes 
        @see axial_flux  
        """
        n = len(self.rs.segments)
        return np.array([self.axial_flux(i, sim_time, rx, sxx, k_soil, cells, True) for i in range(0, n)])
        
    def radial_fluxes(self, sim_time, rx, sxx, k_soil=[], cells=True):
        """ returns the exact radial fluxes (calls base class)
            @param sim_time [day]       needed for age dependent conductivities (age = sim_time - segment creation time)        
            @param rx [cm]              root xylem matric potentials per root system node
            @param sxx [cm]             soil matric potentials given per segment or per soil cell
            @param k_soil [day-1]       optionally, soil conductivities can be prescribed per segment, 
                                        conductivity at the root surface will be limited by the value, i.e. kr = min(kr_root, k_soil) 
            @param cells                indicates if the matric potentials are given per cell (True) or by segments (False)
            @return [cm3 day-1] radial volumetric flow rate            
        """
        return np.array(self.segFluxes(sim_time, rx, sxx, False, cells, k_soil))  # approx = False
        
    def get_nodes(self):
        """ converts the list of Vector3d to a 2D numpy array """
        return np.array(list(map(lambda x: np.array(x), self.rs.nodes)))

    def get_segments(self):
        """ converts the list of Vector2i to a 2D numpy array """
        return np.array(list(map(lambda x: np.array(x), self.rs.segments)), dtype=np.int64)

    def get_ages(self, final_age=0.):
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
        np.sort(rootsegments, axis=None)	
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
        get_y_node = lambda vec : vec[1]
        get_x_node = lambda vec : vec[0]
        get_nodetype = lambda y : organTypes[y - 1]
        nodesy = np.array([get_y_node(xi) for xi in segments], dtype=np.int64)
        nodesx = np.array([get_x_node(xi) for xi in segments], dtype=np.int64)
        nodesy = np.setdiff1d(nodesy, nodesx)  # select all the nodes which belong to tip of an organ
        nodes_type = np.array([get_nodetype(xi) for xi in nodesy], dtype=np.int64)
        tiproots = np.intersect1d(np.where(nodes_type == 2, nodesy, -1), nodesy)  # take root tips
        tipstem = np.intersect1d(np.where(nodes_type == 3, nodesy, -1), nodesy)  # take stem tips
        tipleaf = np.intersect1d(np.where(nodes_type == 4, nodesy, -1), nodesy)  # take leaf tips
        return tiproots, tipstem, tipleaf 
        
    def get_organ_segments_tips(self):	
        """ return index of segments at the end of each organ """
        tiproots, tipstems, tipleaves = self.get_organ_nodes_tips() 
        tiproots = tiproots - np.ones(tiproots.shape, dtype=np.int64)  # segIndx = seg.y -1
        tipstems = tipstems - np.ones(tipstems.shape, dtype=np.int64)  # segIndx = seg.y -1
        tipleaves = tipleaves - np.ones(tipleaves.shape, dtype=np.int64)  # segIndx = seg.y -1
        return tiproots, tipstems, tipleaves
        
    def test(self):
        """ perfoms some sanity checks, and prints to the console """
        print("\nXylemFluxPython.test:")
        # 1 check if segment index is node index-1
        segments = self.get_segments()
        for i, s_ in enumerate(segments):
            if i != s_[1] - 1:
                raise "Segment indices are mixed up"
        print(len(segments), "segments")
        # 2 check for very small segments
        seg_length = self.segLength()
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
        print("{:g} different root types".format(np.max(types) + 1))
        # 4 Print segment age range
        ages = self.get_ages()
        print("ages from {:g} to {:g}".format(np.min(ages), np.max(ages)))
        print()
        
    def summarize_fluxes(self, fluxes, sim_time, rx, p_s, k_soil=[], cells=False, show_matrices=False):
        """gives an overview of the radial and axial water flux. allows us to check that the water balance is about 0
            @param sim_time [day]       needed for age dependent conductivities (age = sim_time - segment creation time)        
            @param rx [cm]              root xylem matric potentials per root system node
            @param p_s [cm]             outer matric potentials given per segment or per cell
            @param k_soil [day-1]       optionally, soil conductivities can be prescribed per segment, 
                                        conductivity at the root surface will be limited by the value, i.e. kr = min(kr_root, k_soil) 
            @param cells                indicates if the matric potentials are given per cell (True) or by segments (False)
            @param show_matrices        show water flux matrices
            @return [cm3 day-1] radial volumetric flow rate            
        
        """
        organTypes = np.array(self.rs.organTypes)
        sumfluxes = [sum(np.where(organTypes == 2, fluxes, 0)), sum(np.where(organTypes == 3, fluxes, 0)), sum(np.where(organTypes == 4, fluxes, 0))]	
        tipsegroots, tipsegstems, tipsegleaves = self.get_organ_segments_tips() 
        axialfluxes_j = np.array([self.axial_flux(i, sim_time, rx=rx, sxx=p_s, k_soil=k_soil, ij=True, cells=cells) for i in range(0, len(self.rs.segments))]) 
        axialfluxes_i = -np.array([self.axial_flux(i, sim_time, rx=rx, sxx=p_s, k_soil=k_soil, ij=False, cells=cells) for i in range(0, len(self.rs.segments))])
        
        balance = -axialfluxes_i - axialfluxes_j + fluxes
        if show_matrices :
            print("matrix of radial fluxes per segments  (<0 leaving segment) :")
            print(-np.array(fluxes))
            print("axial water flux at segment i node (<0 leaving segment) :")
            print(axialfluxes_j)
            print("axial water flux at segment j node (<0 leaving segment) :")
            print(axialfluxes_i)
            print("water balance for each segment (sum of axial and radial water fluxes): ")
            print(balance)
        
        print("\nradial water fluxes (<0 leaving segments)")
        print('radial water flux for roots is {} [cm3/day]'.format(-sumfluxes[0]))
        print('radial water flux for stems is {} [cm3/day]'.format(-sumfluxes[1]))
        print('radial water flux for leaves is {} [cm3/day]'.format(-sumfluxes[2]))
        print('net radial water flux is {} [cm3/day]'.format(-sum(fluxes)))
        print("\naxial water flux at the tip of organs (<0 leaving segments)")
        print('axial water flux at the tip of the roots is {} [cm3/day]'.format(-sum(axialfluxes_i[tipsegroots])))
        if len(tipsegstems) > 0:
            print('axial water flux at the tip of the stems is {} [cm3/day]'.format(-sum(axialfluxes_j[tipsegstems])))
        if len(tipsegleaves) > 0:
            print('axial water flux at the tip of the leaves is {} [cm3/day]'.format(-sum(axialfluxes_j[tipsegleaves])))
        print('sum (absolut value) of net water balance in each segment {} [cm3/day] (should be 0)'.format(sum(abs(balance)))) 
        
        # check that for each node, influx = out flux:
        organTypes = self.get_organ_types()
        segments = self.get_segments() 	
        get_y_node = lambda vec : np.array(vec)[1]
        get_x_node = lambda vec : np.array(vec)[0]
        get_nodetype = lambda y : organTypes[y - 1]
        nodesy = np.array([get_y_node(xi) for xi in segments])
        nodesx = np.array([get_x_node(xi) for xi in segments])
        nodeflux = np.full(len(self.get_nodes()), math.nan)
        nodes = np.arange(len(self.get_nodes()))
        for nd in range(len(nodeflux)) :
            fluxseg_roots_bellow = axialfluxes_j[np.logical_and((nodesx == nd), (organTypes == 2))]  # upper root flux where node is upper node
            fluxseg_roots_above = axialfluxes_i[np.logical_and((nodesy == nd), (organTypes == 2))]  # lower root flux where node is lower node
            fluxseg_stemsleaves_bellow = axialfluxes_j[np.logical_and((nodesy == nd), (organTypes > 2))]  # upper root flux where node is upper node
            fluxseg_stemsleaves_above = axialfluxes_i[np.logical_and((nodesx == nd), (organTypes > 2))]  # lower root flux where node is lower node
            nodeflux[nd] = sum((sum(fluxseg_roots_bellow), sum(fluxseg_roots_above), sum(fluxseg_stemsleaves_bellow), sum(fluxseg_stemsleaves_above)))
        tiproots, tipstems, tipleaves = self.get_organ_nodes_tips() 
        tips = np.concatenate((tiproots, tipstems, tipleaves))
        nodeflux[tips] = math.nan
        print('sum (absolut value) of net water balance in each node, appart from node at the tip of organs {} [cm3/day] (should be 0)'.format(np.nansum(abs(nodeflux))))
        tot_balance = -sum(axialfluxes_i[tipsegroots]) - sum(axialfluxes_j[tipsegstems]) - sum(axialfluxes_j[tipsegleaves]) - sum(fluxes)
        print('net water flux in plant: {} [cm3/day] (should be 0)'.format(tot_balance))
        return fluxes

    @staticmethod
    def read_rsml(file_name :str):
        """ reads an RSML file and converts to MappedSegments with units [cm]
        @file_name     the file name of the rsml, including file extension (e.g. "test.rsml" ) 
        @return a CPlantBox MappedSegments object
        """
        polylines, props, funcs = rsml.read_rsml(file_name)
        bn = 0  # count base roots
        for i, _ in enumerate(polylines):
            if props["parent-poly"][i] < 0:
                bn += 1
        if bn > 1: 
            print("XylemFlux.read_rsml: added an artificial shoot")
            polylines, props, funcs = rsml.artificial_shoot(polylines, props, funcs)        
        nodes, segs = rsml.get_segments(polylines, props)
        radii, seg_ct, types = rsml.get_parameter(polylines, funcs, props)
        print("XylemFluxPython.read_rsml: read rsml with", len(nodes), "nodes and", len(radii), "radii")
        nodes = np.array(nodes)  # for slicing in the plots
        nodes2 = []  # Conversions...
        for n in nodes:
            nodes2.append(pb.Vector3d(n[0] , n[1] , n[2]))
        segs2 = []
        nodeCTs = np.zeros((len(nodes), 1))  # we need node creation times
        for i, s in enumerate(segs):
            nodeCTs[s[1]] = seg_ct[i]
            segs2.append(pb.Vector2i(int(s[0]), int(s[1])))
        radii = np.array(radii)
        types = np.array(types, dtype=np.int64) - 1  # index must start with 0
        return pb.MappedSegments(nodes2, nodeCTs, segs2, radii, types)  # root system grid

    @staticmethod
    def zero_rows(M, rows):
        """ Sets the given rows of the sparse matrix M to zero (I do not know the best way to do it)
        @param M         sparse matrix
        @param rows      list of row indices
        """
        diag = sparse.eye(M.shape[0]).tolil()
        for r in rows:
            diag[r, r] = 0
        return diag.dot(M)

    @staticmethod
    def bc_dirichlet(Q, b, n0, d):
        """ prescribes a Dirichlet boundary conditions for the system Qx=b
        @param Q          system matrix
        @param b          rhs vector
        @param n0         list of node indices, where the Dirichlet bc is applied
        @param d [cm]     list of Dirichlet values   
        @return Q, b, the updated matrix, and rhs vector 
        """
        assert len(n0) == len(d), "XylemFlux.bc_dirichlet: number of nodes n0 and Dirichlet values d must be equal"
        Q = XylemFluxPython.zero_rows(Q, n0)
        for c in range(0, len(n0)):
            i = n0[c]
            Q[i, i] = 1
            b[i] = d[c]
        return Q, b

    @staticmethod
    def bc_neumann(Q, b, n0, f):
        """ prescribes a Neumann boundary condition for the root system Qx=b
        @param Q                      system matrix
        @param b                      rhs vector
        @param n0                     list of node indices, where the Neumann boundary condition is applied
        @param f [cm3 day-1]          list of Neumann values   
        @return Q, b, the updated matrix, and rhs vector                 
        """
        for c in range(0, len(n0)):
            b[int(n0[c])] += f[c]
        return Q, b

    @staticmethod
    def convert_(x, dtype=np.float64):
        """ not used anymore (?) """
        return np.array(list(map(lambda x: np.array(x, dtype), x)), dtype)  # is there a better way?
    
