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
        
        self.seg_ind = [0] # segment indices for Neuman flux
        self.node_ind = [0] # node indices for Dirichlet flux

    def solve_neumann(self, sim_time :float, value, sxx, cells :bool, soil_k = []) :
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
        if isinstance(value, (float,int)):
            n = len(self.seg_ind)
            value = [value/n]*n
            
        if len(soil_k) > 0:
            self.linearSystem(sim_time, sxx, cells, soil_k) # C++ (see XylemFlux.cpp)
        else:
            self.linearSystem(sim_time, sxx, cells) # C++ (see XylemFlux.cpp)
        Q = sparse.coo_matrix((np.array(self.aV), (np.array(self.aI), np.array(self.aJ))))
        Q = sparse.csr_matrix(Q)
        Q, b = self.bc_neumann(Q, self.aB, self.seg_ind, value)  # cm3 day-1
        x = LA.spsolve(Q, b, use_umfpack = True)  # direct
        # print ("linear system assembled and solved in", timeit.default_timer() - start, " s")
        return x

    def solve_dirichlet(self, sim_time :float, value :list, sxc :float, sxx, cells :bool, soil_k = []):
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
        if isinstance(value, (float,int)):
            n = len(self.node_ind)
            value = [value]*n       
         
        if len(soil_k) > 0:
            self.linearSystem(sim_time, sxx, cells, soil_k)  # C++ (see XylemFlux.cpp)
        else:
            self.linearSystem(sim_time, sxx, cells)
        self.linearSystem(sim_time, sxx, cells)  # C++ (see XylemFlux.cpp)
        Q = sparse.coo_matrix((np.array(self.aV), (np.array(self.aI), np.array(self.aJ))))
        Q = sparse.csr_matrix(Q)
        Q, b = self.bc_dirichlet(Q, self.aB, self.node_ind, value)
        x = LA.spsolve(Q, b, use_umfpack = True)
        return x

    def solve(self, sim_time :float, trans :list, sx :float, sxx, cells :bool, wilting_point :float, soil_k = []):
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

            x = self.solve_neumann(sim_time, trans, sxx, cells, soil_k) # try neumann, if below wilting point, switch to Dirichlet

            if x[0] <= wilting_point:
                Q = sparse.coo_matrix((np.array(self.aV), (np.array(self.aI), np.array(self.aJ))))
                Q = sparse.csr_matrix(Q)
                Q, b = self.bc_dirichlet(Q, self.aB, [0], [float(wilting_point)])
                try:
                    x = LA.spsolve(Q, b, use_umfpack = True)
                except:
                    print("XylemFluxPython.solve: Exeption solving Dirichlet")
                    print("Dirichlet at ", trans - sx, "cm")
                    print("b", b)
                    raise
        else:
            print("XylemFluxPython.solve: used Dirichlet because collar cell soil matric potential is below wilting point", sx)
            x = self.solve_dirichlet(sim_time, wilting_point, sx, sxx, cells, soil_k)

        return x
    
    
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
        s = self.rs.segments[seg_ind]  # collar segment
        if len(k_soil) > 0:
            ksoil = k_soil[seg_ind]
        else:
            ksoil = 1.e9 # much
        if ij: 
            i, j = int(s.x), int(s.y)  # node indices
        else:
            j, i = int(s.x), int(s.y)  
        n1, n2 = self.rs.nodes[i], self.rs.nodes[j]  # nodes
        v = n2.minus(n1)
        l = v.length()  # length of segment
        v.normalize()  # normalized v.z is needed for qz
        if cells:
            cell_ind = self.rs.seg2cell[seg_ind]
            if cell_ind>=0:
                p_s = sxx[cell_ind]  # soil pressure at collar segment
            else:
                p_s = self.airPressure
        else:
            p_s = sxx[seg_ind]
        a = self.rs.radii[seg_ind]  # radius
        ot = int(self.rs.organTypes[seg_ind])  # conductivities kr, kx
        st = int(self.rs.subTypes[seg_ind])  # conductivities kr, kx
        age = sim_time - self.rs.nodeCTs[int(s.y)]
        kr = self.kr_f(age, st, ot)  # c++ conductivity call back functions
        kr = min(kr, ksoil)
        kx = self.kx_f(age, st, ot)
        tau = math.sqrt(2 * a * math.pi * kr / kx)  # cm-2
        AA = np.array([[1, 1], [math.exp(tau* l), math.exp(-tau* l)] ])          
        bb = np.array([rx[i] - p_s, rx[j] - p_s])  # solve for solution
        d = np.linalg.solve(AA, bb)  # compute constants d_1 and d_2 from bc
        dpdz0 = d[0] * tau - d[1] * tau # insert z = 0, z = l into exact solution   
        f = kx * (dpdz0 + v.z)
        if ij:
            f = f*(-1)
        return f  

    def collar_flux(self, sim_time, rx, sxx, k_soil = [], cells = True):
        """ returns the exact transpirational flux of the xylem model solution @param rx
            @see axial_flux        
        """
        return self.axial_flux(0, sim_time, rx, sxx, k_soil, cells, True)

    def axial_fluxes(self, sim_time, rx, sxx, k_soil = [], cells = True):
        """ returns the axial fluxes 
        @see axial_flux  
        """
        n = len(self.rs.segments)
        return np.array([self.axial_flux(i, sim_time, rx, sxx, k_soil, cells, True) for i in range(0, n)])
        
    def radial_fluxes(self, sim_time, rx, sxx, k_soil = [], cells = True):
        """ returns the exact radial fluxes (calls base class)
            @param sim_time [day]       needed for age dependent conductivities (age = sim_time - segment creation time)        
            @param rx [cm]              root xylem matric potentials per root system node
            @param sxx [cm]             soil matric potentials given per segment or per soil cell
            @param k_soil [day-1]       optionally, soil conductivities can be prescribed per segment, 
                                        conductivity at the root surface will be limited by the value, i.e. kr = min(kr_root, k_soil) 
            @param cells                indicates if the matric potentials are given per cell (True) or by segments (False)
            @return [cm3 day-1] radial volumetric flow rate            
        """
        return np.array(self.segFluxes(sim_time, rx, sxx, False, cells, k_soil)) # approx = False
        
    def get_outer_matpot_matix(self, p_s, p_a):
        # do we stil need that one? could be replaced by return p_s*(ot==2) + p_a*(ot!=2) in the script
        p_out = self.get_organ_types()	
        for i in range(0, len(p_out)):
            if(p_out[i] == 2): #= root type
                p_out[i] = p_s
            else: #= stem or leaf type
                p_out[i] = p_a
        return p_out
        
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
        
    def get_nodes_index(self,ot):	
        """ return node indices of segments with organ type @param ot """
        segments = self.get_segments() 	
        nodes = self.get_nodes()	
        organTypes = self.get_organ_types()	
        rootsegments = segments[organTypes == ot]	
        rootsegments.flatten()	
        np.sort(rootsegments, axis=None)	
        nodesidx =  np.unique(rootsegments)	
        return nodesidx	
        	
    def get_nodes_organ_type(self,ot):	
        """ return node coordinates of segments with organ type @param ot """
#         nodes = self.get_nodes()
#         return nodes[self.get_nodes_index()] # should do the job?       
        segments = self.get_segments() 	
        nodes = self.get_nodes()	
        organTypes = self.get_organ_types()	
        rootsegments = segments[organTypes == ot]	
        rootsegments.flatten()	
        np.sort(rootsegments, axis=None)	
        rootnodes =  np.unique(rootsegments)	
        nodes = nodes[rootnodes]	
        return nodes
        
    def get_segments_index(self,ot):	
        """ return node indices of organ type @param ot """
        organTypes = self.get_organ_types()
        segIdx = np.array( list(range(0, len(organTypes))))
        otsegs = segIdx[organTypes == ot]	
        return otsegs
        
    def get_organ_types(self):
        """ segment organ types as numpy array """
        return np.array(self.rs.organTypes)
        
    def get_subtypes(self):
        """ segment sub types as numpy array """
        return np.array(self.rs.subTypes)
        
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
            raise "Types start with an index greater than 0"
        print("{:g} different root types".format(np.max(types) + 1))
        # 4 Print segment age range
        ages = self.get_ages()
        print("ages from {:g} to {:g}".format(np.min(ages), np.max(ages)))
        print()

    @staticmethod
    def read_rsml(file_name :str):
        """ reads an RSML file and converts to MappedSegments with units [cm]
        @file_name     the file name of the rsml, including file extension (e.g. "test.rsml" ) 
        @return a CPlantBox MappedSegments object
        """
        polylines, props, funcs = rsml.read_rsml(file_name)
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
        types = np.array(types, dtype = np.int64) - 1  # index must start with 0
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
        assert len(n0)==len(d), "XylemFlux.bc_dirichlet: number of nodes n0 and Dirichlet values d must be equal"
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
            b[n0[c]] += f[c]
        return Q, b

    @staticmethod
    def convert_(x, dtype = np.float64):
        """ not used anymore (?) """
        return np.array(list(map(lambda x: np.array(x, dtype), x)), dtype)  # is there a better way?
    