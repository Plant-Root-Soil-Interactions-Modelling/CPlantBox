import timeit

import types
import numpy as np
from scipy import sparse
import scipy.sparse.linalg as LA
import matplotlib.pyplot as plt

import plantbox as pb
from plantbox import PlantHydraulicModel as PlantHydraulicModelCPP
from structural.MappedOrganism import MappedPlantPython
from functional.Perirhizal import PerirhizalPython as Perirhizal

import rsml.rsml_reader as rsml


class PlantHydraulicModel(PlantHydraulicModelCPP):
    """  Plant hydraulic model super class 
    
        Currently implentations are located in the specialisations below:
            HydraulicModel_Meunier: Meunier hybrid solver (Meunier et al. 2017)    
            HydraulicModel_Doussan: Doussan solver (Doussan et al. 2006)    
            
        To implement a hydraulic model, overwrite: 
            solve_dirichlet
            solve_neumann
            solve
            solve_again (optionally)
            radial fluxes
            axial fluxes
            
        all function use matric potentials as units for input and output 
        
        if the linear system does not change for several solve calls, 
        it is possible to use cached factorizations, by stating cached = True in the constructor
        
        Preliminary tests in 
        dumux-rosi/python/roots/xylem_m31_new.py, xylem_m32_new.py
        dumux-rosi/python/coupled/C12new/coupled_c12_sra.py
    """

    def __init__(self, ms, params, cached = False):
        """ 
        @param ms is of type MappedSegments (or specializations), or a string containing a rsml filename
        @param params hydraulic conductivities described by PlantHydraulicParameters
        @param if cached == "True" sparse factorization is cached for faster solving using solve_again() 
        """
        if isinstance(ms, str):
            ms = self.read_rsml(ms)
            super().__init__(ms, params)
        else:
            super().__init__(ms, params)

        self.ot_root = int(pb.OrganTypes.root)
        self.cached = cached
        self.last = "none"  # after first solve() call "neumann" or "dirichlet"
        self.neumann_ind = [0]  # node indices for Neumann flux
        self.dirichlet_ind = [0]  # node indices for Dirichlet flux
        self.collar_index_ = self.collar_index()  # segment index of the collar segement
        self.wilting_point = -15000  # [cm]

    def solve_dirichlet(self, sim_time:float, collar_pot:list, rsx, cells:bool):
        """ solves the flux equations, with a neumann boundary condtion, see solve()            
            @param sim_time [day]           needed for age dependent conductivities (age = sim_time - segment creation time)
            @param collar_pot [cm3 day-1]   collar potential
            @param rsx [cm]                 soil matric potentials given per segment or per soil cell            
            @param cells                    indicates if the matric potentials are given per cell (True) or by segments (False)              
            @return [cm] root matric potential per root system node         
        """
        raise "PlantHydraulicModel(): use implementations of this abstract super class, e.g. HydraulicModel_Meunier(), or HydraulicModel_Meunier()"

    def solve_neumann(self, sim_time:float, t_act:list, rsx, cells:bool):
        """ solves the flux equations, with a neumann boundary condtion, see solve()
            @param sim_time [day]       needed for age dependent conductivities (age = sim_time - segment creation time)
            @param t_act [cm3 day-1]    tranpirational flux is negative
            @param rsx [cm]             soil matric potentials given per segment or per soil cell
            @param cells                indicates if the matric potentials are given per cell (True) or by segments (False)              
            @return [cm] root matric potential per root system node          
        """
        raise "PlantHydraulicModel(): use implementations of this abstract super class, e.g. HydraulicModel_Meunier(), or HydraulicModel_Meunier()"

    def solve(self, sim_time:float, t_act:list, rsx, cells:bool):
        """ Solves the hydraulic model using Neumann boundary conditions and switching to Dirichlet in case wilting point is reached
            @param sim_time [day]        needed for age dependent conductivities (age = sim_time - segment creation time)
            @param t_act [cm3 day-1]     transpiration rate
            @param rsx [cm]              soil matric potentials given per segment or per soil cell
            @param cells                 indicates if the matric potentials are given per cell (True) or by segments (False)
            @return [cm] root matric potential per root system node  
        """
        raise "PlantHydraulicModel(): use implementations of this abstract super class, e.g. HydraulicModel_Meunier(), or HydraulicModel_Meunier()"

    def solve_again(self, sim_time:float, t_act:list, rsx, cells:bool, soil_k = []):
        """ Solves the hydraulic model using Neumann boundary conditions and switching to Dirichlet in case wilting point is reached
            Depending of the method solve_again() is much faster using chached factorization of the last solve() command
            @param sim_time [day]        needed for age dependent conductivities (age = sim_time - segment creation time)
            @param t_act [cm3 day-1]     transpiration rate
            @param rsx [cm]              soil matric potentials given per segment or per soil cell
            @param cells                 indicates if the matric potentials are given per cell (True) or by segments (False)
            @return [cm] root matric potential per root system node  
        """
        self.solve(sim_time, t_act, rsx, cells, soil_k)

    def radial_fluxes(self, sim_time:float, rx, rsx, cells = False):
        """ returns the radial fluxes per segment [cm3 day-1]"""
        raise "PlantHydraulicModel(): use implementations of this abstract super class, e.g. HydraulicModel_Meunier(), or HydraulicModel_Meunier()"

    def axial_fluxes(self, sim_time:float, rx, rsx, cells = False):
        """ returns the axial fluxes per segment [cm3 day-1]"""
        raise "PlantHydraulicModel(): use implementations of this abstract super class, e.g. HydraulicModel_Meunier(), or HydraulicModel_Meunier()"

    def soil_fluxes(self, sim_time:float, rx, rsx):
        """ sums the radial fluxes over the soil cells, returns a dictionary with sources and sinks (with cell id as key) """
        fluxes = self.radial_fluxes(sim_time, rx, rsx, True)
        return self.sumSegFluxes(fluxes)

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
            print("                           subTypes Doussan[{:g}, {:g}] ".format(np.min(types), np.max(types)))
        return pb.MappedSegments(nodes2, cts, segs2, segRadii, segTypes)  # root system grid

    def collar_index(self):
        """ returns the segment index of the collar segment, node index of the collar node is always 0 """
        segs = self.ms.segments
        for i, s in enumerate(segs):
            if s.x == 0:
                return i

    @staticmethod
    def sinusoidal(t):
        """ sinusoidal function (used for transpiration, integral over one day is 1)"""
        return np.sin(2. * np.pi * np.array(t) - 0.5 * np.pi) + 1.

    @staticmethod
    def sinusoidal2(t, dt):
        """ sinusoidal function from 6:00 - 18:00, 0 otherwise (used for transpiration, integral over one day is 1)"""
        return np.maximum(0., np.pi * (np.cos(2 * np.pi * (t - 0.5)) + np.cos(2 * np.pi * ((t + dt) - 0.5))) / 2)

    def get_nodes(self):
        """ converts the list of Vector3d to a 2D numpy array (from MappedOrganism) """
        nodes = self.ms.nodes
        return np.array(list(map(lambda x: np.array(x), nodes)))

    def get_segments(self):
        """ converts the list of Vector2i to a 2D numpy array """
        segments = self.ms.segments
        return np.array(list(map(lambda x: np.array(x), segments)), dtype = np.int64)

    def get_organ_types(self):
        """ segment organ types as numpy array """
        return np.array(self.ms.organTypes)

    def get_ages(self, final_age = -1.):
        """ ages per segment
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

    def get_efffective_kr(self, sim_time):
        """ effective radial conductivities per segment (radial conductivities multiplied by segment surface) [cm2 day-1] """
        return np.array(self.params.getEff(sim_time))  #

    def get_kr(self, sim_time):
        """ radial conductivities per segment [1 day-1] """
        return np.array(self.params.getKr(sim_time))

    def get_kx(self, sim_time):
        """ axial conductivities per segment [cm3 day-1]"""
        return np.array(self.params.getKx(sim_time))

    def get_hs(self, sx):
        """ soil matric potential per segment [cm] """
        return np.array(self.ms.getHs(sx))

    def get_transpiration(self, sim_time, rx, rsx, cells = False):
        """ actual transpiration [cm3 day-1], calculated as the sum of all radial fluxes"""
        return np.sum(self.radial_fluxes(sim_time, rx, rsx, cells))

    def get_krs(self, sim_time):
        """ calculatets root system conductivity [cm2/day] at simulation time @param sim_time [day]
        if there is no single collar segment at index 0, pass indices using @param seg_ind, see find_base_segments
        """
        n = self.ms.getNumberOfMappedSegments()
        rsx = np.ones((n, 1)) * (-500)
        rsx = self.ms.total2matric(rsx)
        rx = self.solve_dirichlet(sim_time, -15000, rsx, cells = False)
        t_act = -self.get_transpiration(sim_time, rx, rsx, cells = False)
        nodes = self.ms.nodes  # bit heavy
        segs = self.ms.segments
        krs = t_act / (-500 - 0.5 * (nodes[segs[0].x].z + nodes[segs[0].y].z) - rx[0])
        return krs , t_act

    def get_suf(self, sim_time):
        """ Standard uptake fraction (SUF) [1] per root segment, should add up to 1 """
        n = self.ms.getNumberOfMappedSegments()
        rsx = np.ones((n, 1)) * (-500)
        rsx = self.ms.total2matric(rsx)
        rx = self.solve_dirichlet(sim_time, -15000, rsx, cells = False)
        q = self.radial_fluxes(sim_time, rx, rsx)
        return np.array(q) / np.sum(q)

    def get_heff(self, rsx):
        """ effective total potential [cm] """
        suf = self.get_suf()
        heff = suf.dot(self.ms.matric2total(rsx))
        return heff[0]

    def test_unmapped(self):
        """ test for unmapped segments """
        pass_ = True
        segs = self.get_segments()
        types = self.ms.subTypes
        ots = self.ms.organTypes
        nodes = self.get_nodes()
        map = self.ms.seg2cell
        for seg_id, cell_id in map.items():
            if cell_id < 0:
                print("Warning: segment ", seg_id, "is not mapped, this will cause problems with coupling!",
                      "organType", ots[seg_id], "subType", types[seg_id], "mid", nodes[segs[seg_id][0]], nodes[segs[seg_id][1]])
                pass_ = False
        return pass_

    def test(self):
        """ perfoms some sanity checks, and prints to the console """
        print("PlantHydraulicModel.test():")
        # 1 check if segment index is node index-1
        segments = self.get_segments()
        nodes = self.get_nodes()
        types = self.ms.subTypes
        for i, s_ in enumerate(segments):
            if i != s_[1] - 1:
                raise "Error: Segment indices are mixed up!"
        print(len(nodes), "nodes:")
        for i in range(0, min(5, len(nodes))):
            print("Node", i, nodes[i])
        print(len(segments), "segments:")
        # 1b check if there are multiple basal roots (TODO)
        for i in range(0, min(5, len(segments))):
            print("Segment", i, segments[i], "subType", types[i])
        ci = self.collar_index()
        self.collar_index_ = ci
        print("Collar segment index", ci)
        print("Collar segment", segments[ci])
        first = True
        for s in segments:
            if s[0] == 0:
                if first:
                    first = False
                else:
                    print("Warning: multiple segments emerge from collar node (always node index 0)", ci, s)
        # 2 check for very small segments
        seg_length = self.ms.segLength()
        c = 0
        for i, l in enumerate(seg_length):
            if l < 1e-5:
                print(i, l)
                c += 1
        print(c, "segments with length < 1.e-5 cm")
        # 3 check for type range, index should start at 0
        if np.min(types) > 0:
            print("Warning: types start with index", np.min(types), "> 0 !")
        print("{:g} different root types from {:g} to {:g}".format(np.max(types) - np.min(types) + 1, np.min(types), np.max(types)))
        # 4 Print segment age range
        ages = self.get_ages()
        print("ages from {:g} to {:g}".format(np.min(ages), np.max(ages)))
        # 4 check for unmapped indices
        print("segments", len(segments), len(self.ms.segments))
        map = self.ms.seg2cell
        for seg_id, cell_id in map.items():
            # if seg_id < 5:
            #     print("seg_id", seg_id, "cell_id", cell_id)
            if (cell_id < 0) & (types[seg_id] == int(pb.root)):
                print("Warning: root segment ", seg_id, "is not mapped (out of the soil domain), this could cause problems with coupling!",
                nodes[segments[seg_id][0]], nodes[segments[seg_id][1]])
        print()


class HydraulicModel_Meunier(PlantHydraulicModel):
    """
    Meunier hybrid solver (Meunier et al. 2017) 
    
    For perfomance reasons building the matrix is implented in C++ (PlantHydraulicModelCPP)
    
    If cached is true, for each solve() call a sparse LU factorization is stored, which is used in every solve_again() call
    (no convenient sparse cholesky implementation in scipy)
    """

    def __init__(self, ms, params, cached = False):
        """ 
            @param ms is of type MappedSegments (or specializations), or a string containing a rsml filename
            @param params hydraulic conductivities described by PlantHydraulicParameters
            @param if cached == "True" sparse factorization is cached for faster solving using solve_again() 
        """
        super().__init__(ms, params, cached)
        self.usecached_ = False

    def solve_dirichlet(self, sim_time:float, collar_pot:list, rsx, cells:bool, soil_k = []):
        """ solves the flux equations, with a neumann boundary condtion, see solve()
            @param sim_time [day]           needed for age dependent conductivities (age = sim_time - segment creation time)
            @param collar_pot [cm]          collar potential
            @param rsx [cm]                 soil matric potentials given per segment or per soil cell            
            @param cells                    indicates if the matric potentials are given per cell (True) or by segments (False) 
            @param soil_k [day-1]     optionally, soil conductivities can be prescribed per segment, 
                                      conductivity at the root surface will be limited by the value, i.e. kr = min(kr_root, k_soil)   
            @return [cm] root matric potential per root system node         
        """
        self.last = "dirichlet"
        self.linearSystemMeunier(sim_time, rsx, cells, soil_k)
        Q = sparse.csc_matrix((np.array(self.aV), (np.array(self.aI), np.array(self.aJ))))

        if isinstance(collar_pot, (float, int)):
            collar_pot = [collar_pot] * len(self.dirichlet_ind)

        Q, b = self.bc_dirichlet(Q, self.aB, self.dirichlet_ind, collar_pot)

        if self.usecached_:  # TODO only self.aB is needed;  building rhs b could be improved
            return self.dirichletB.solve(np.array(b))
        else:
            return LA.spsolve(Q, b, use_umfpack = True)

    def solve_neumann(self, sim_time:float, t_act:list, rsx, cells:bool, soil_k = []):
        """ solves the flux equations, with a neumann boundary condtion, see solve()
            @param sim_time [day]       needed for age dependent conductivities (age = sim_time - segment creation time)
            @param t_act [cm3 day-1]    tranpirational flux is negative
            @param rsx [cm]             soil matric potentials given per segment or per soil cell
            @param cells                indicates if the matric potentials are given per cell (True) or by segments (False)  
            @param soil_k [day-1]     optionally, soil conductivities can be prescribed per segment, 
                                      conductivity at the root surface will be limited by the value, i.e. kr = min(kr_root, k_soil)  
            @return [cm] root matric potential per root system node         
        """
        self.last = "neumann"
        self.linearSystemMeunier(sim_time, rsx, cells, soil_k)  # C++ (see XylemFlux.cpp)
        Q = sparse.csc_matrix((np.array(self.aV), (np.array(self.aI), np.array(self.aJ))))

        if isinstance(t_act, (float, int)):
            t_act = [t_act / len(self.neumann_ind)] * len(self.neumann_ind)

        b = self.bc_neumann(self.aB, self.neumann_ind, t_act)

        if self.usecached_:
            return self.neumannB.solve(np.array(b))
        else:
            return LA.spsolve(Q, b, use_umfpack = True)

    def solve(self, sim_time:float, t_act:list, rsx, cells:bool, soil_k = []):
        """ Solves the hydraulic model using Neumann boundary conditions and switching to Dirichlet in case wilting point is reached
            @param sim_time [day]        needed for age dependent conductivities (age = sim_time - segment creation time)
            @param t_act [cm3 day-1]     transpiration rate
            @param rsx [cm]              soil matric potentials given per segment or per soil cell
            @param cells                 indicates if the matric potentials are given per cell (True) or by segments (False)
            @param soil_k [day-1]     optionally, soil conductivities can be prescribed per segment, 
                                      conductivity at the root surface will be limited by the value, i.e. kr = min(kr_root, k_soil)  
            @return [cm] root matric potential per root system node  
        """
        if self.cached:  # store sparse LU factorization
            self.linearSystemMeunier(sim_time, rsx, cells, soil_k)  # self.aV, self.aB, self.aI, self.aJ
            Q = sparse.csc_matrix(sparse.coo_matrix((np.array(self.aV), (np.array(self.aI), np.array(self.aJ)))))
            self.neumannB = LA.splu(Q)  # for Neumann Q
            Q, b = self.bc_dirichlet(Q, self.aB, self.dirichlet_ind, [self.wilting_point] * len(self.dirichlet_ind))
            self.dirichletB = LA.splu(Q)  # for Dirichlet Q
            x = self.solve_again(sim_time, t_act, rsx, cells)
        else:
            x = self.solve_neumann(sim_time, t_act, rsx, cells, soil_k)  # try neumann, if below wilting point, switch to Dirichlet
            self.last = "neumann"
            if x[0] <= self.wilting_point:
                Q = sparse.csc_matrix((np.array(self.aV), (np.array(self.aI), np.array(self.aJ))))
                Q, b = self.bc_dirichlet(Q, self.aB, self.dirichlet_ind, [self.wilting_point] * len(self.dirichlet_ind))
                x = LA.spsolve(Q, b, use_umfpack = True)
                self.last = "dirichlet"

        return x

    def solve_again(self, sim_time:float, t_act:list, rsx, cells:bool, soil_k = []):
        """ Solves the hydraulic model using Neumann boundary conditions and switching to Dirichlet in case wilting point is reached
            Depending of the method solve_again() is much faster using chached factorization of the last solve() command
            @param sim_time [day]        needed for age dependent conductivities (age = sim_time - segment creation time)
            @param t_act [cm3 day-1]     transpiration rate
            @param rsx [cm]              soil matric potentials given per segment or per soil cell
            @param cells                 indicates if the matric potentials are given per cell (True) or by segments (False)
            @param soil_k [day-1]     optionally, soil conductivities can be prescribed per segment, 
                                      conductivity at the root surface will be limited by the value, i.e. kr = min(kr_root, k_soil)  
            @return [cm] root matric potential per root system node  
        """
        if not self.cached:
            print("HydraulicModel_Meunier.solve_again() warning: call without cached==True")
            return self.solve(sim_time, t_act, rsx, cells)
        self.usecached_ = True  # makes solve_neumann() and solve_dirichlet() to use the precomputed LU factorizuation
        x = self.solve_neumann(sim_time, t_act, rsx, cells, soil_k)  # try neumann, if below wilting point, switch to Dirichlet
        self.last = "neumann"
        self.usecached_ = False  # only True, during solve_again() call
        if x[0] <= self.wilting_point:
            Q = sparse.csc_matrix((np.array(self.aV), (np.array(self.aI), np.array(self.aJ))))
            Q, b = self.bc_dirichlet(Q, self.aB, self.dirichlet_ind, [self.wilting_point] * len(self.dirichlet_ind))  # TODO improve
            x = self.dirichletB.solve(np.array(b))
            self.last = "dirichlet"

        return x

    def get_transpiration(self, sim_time, rx, rsx, cells = False, soil_k = []):
        """ actual transpiration [cm3 day-1], calculated as the sum of all radial fluxes"""
        return np.sum(self.radial_fluxes(sim_time, rx, rsx, cells, soil_k))
        
    def radial_fluxes(self, sim_time, rx, rsx, cells = False, soil_k = []):
        """ returns the exact radial fluxes per segment (calls base class)
            @param sim_time [day]       needed for age dependent conductivities (age = sim_time - segment creation time)        
            @param rx [cm]              root xylem matric potentials per root system node
            @param rsx [cm]             soil matric potentials given per segment or per soil cell
            @param cells                indicates if the matric potentials are given per cell (True) or by segments (False)
            @param soil_k [day-1]     optionally, soil conductivities can be prescribed per segment, 
                                      conductivity at the root surface will be limited by the value, i.e. kr = min(kr_root, k_soil)  
            @return [cm3 day-1] radial volumetric flow rate            
        """
        return np.array(self.getRadialFluxes(sim_time, rx, rsx, False, cells, soil_k))  # approx = False

    def axial_fluxes(self, sim_time, rx, rsx, cells = False):
        """ returns the axial fluxes 
        @see axial_flux  
        """
        n = self.ms.getNumberOfMappedSegments()
        return np.array([self.axial_flux(i, sim_time, rx, rsx, cells, True) for i in range(0, n)])

    def axial_flux(self, seg_ind, sim_time, rx, rsx, cells = True, ij = True):
        """ returns the exact axial flux of segment ij of xylem model solution @param rx
            @param seg_ind              segment index 
            @param sim_time [day]       needed for age dependent conductivities (age = sim_time - segment creation time)        
            @param rx [cm]              root xylem matric potentials per root system node
            @param rsx [cm]             soil matric potentials given per segment or per soil cell
            @param cells                indicates if the matric potentials are given per cell (True) or by segments (False)
            @param ij                   True: calculate axial flux in node i, False: in node j; note that they are not equal due to radial fluxes 
            @return [cm3 day-1] axial volumetric flow rate             
        """
        s = self.ms.segments[seg_ind]
        nodes = self.ms.nodes
        numleaf = 0
        organTypes = self.get_organ_types()
        ot = int(organTypes[seg_ind])  # for conductivities kr, kx
        # node x and y of stem and leave segments are revers with regards to the nodes x and y of roots.
        if sum(((not ij) , (ot == 4) or (ot == 3))) == 1:  # check if only one of the condition is true
            j, i = int(s.x), int(s.y)  # node indices
        else:
            i, j = int(s.x), int(s.y)
        n1, n2 = self.ms.nodes[i], self.ms.nodes[j]  # nodes
        v = n2.minus(n1)
        l = v.length()  # length of segment
        v.normalize()  # normalized v.z is needed for qz
        if cells:
            cell_ind = self.ms.seg2cell[seg_ind]
            if cell_ind >= 0:  # y node belowground
                if len(rsx) > 1:
                    p_s = rsx[cell_ind]  # soil pressure at collar segment
                else:
                    p_s = rsx[0]
            else:
                p_s = self.airPressure
        else:
            p_s = rsx[seg_ind]
        a = self.ms.getEffectiveRadius(seg_ind)  # radius
        st = int(self.ms.subTypes[seg_ind])  # conductivities kr, kx
        age = sim_time - self.ms.nodeCTs[int(s.y)]
        if ot == 4:  # to know which x-th leaf segment it is, to fetch the right gs value
                indices = [i for i, x in enumerate(organTypes) if x == 4]
                numleaf = indices.index(seg_ind)
                if self.pg[0] != 0:
                    p_s = self.pg[numleaf]
        kr = self.params.kr_f(seg_ind, age, st, ot)  # c++ conductivity call back functions
        kx = self.params.kx_f(seg_ind, age, st, ot)
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
        return f

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
        Q = HydraulicModel_Meunier.zero_rows(Q, n0)
        for c in range(0, len(n0)):
            i = n0[c]
            Q[i, i] = 1
            b[i] = d[c]
        return Q, b

    @staticmethod
    def bc_neumann(b, n0, f):
        """ prescribes a Neumann boundary condition for the root system Qx=b
        @param b                      rhs vector
        @param n0                     list of node indices, where the Neumann boundary condition is applied
        @param f [cm3 day-1]          list of Neumann values   
        @return b                     rhs vector                 
        """
        for c in range(0, len(n0)):
            b[int(n0[c])] += f[c]
        return b


class HydraulicModel_Doussan(PlantHydraulicModel):
    """
    Doussan solver (Doussan et al. 2006)  
    """

    def __init__(self, ms, params, cached = True):
        """ 
            @param ms is of type MappedSegments (or specializations), or a string containing a rsml filename
            @param params hydraulic conductivities described by PlantHydraulicParameters
            @param if cached == "True" sparse factorization is cached for faster solving using solve_again() 
        """
        super().__init__(ms, params, cached)
        self.usecached_ = False

    def solve_dirichlet(self, sim_time:float, collar_pot:list, rsx, cells:bool):
        """ solves the flux equations, with a neumann boundary condtion, see solve()            
            @param sim_time [day]           needed for age dependent conductivities (age = sim_time - segment creation time)
            @param collar_pot [cm3 day-1]   collar potential
            @param rsx [cm]                 soil matric potentials given per segment or per soil cell            
            @param cells                    indicates if the matric potentials are given per cell (True) or by segments (False)              
            @return [cm] root matric potential per root system node         
        """
        self.last = "dirichlet"
        if cells:
            rsx = self.get_hs(rsx)  # matric potential per root segment
        if not self.usecached_:
            self.update(sim_time)
        rsx_ = self.ms.matric2total(rsx)
        b = self.Kr.dot(rsx_)
        b[self.collar_index_] += self.kx0 * collar_pot
        rx_ = self.A_d_splu.solve(b)
        rx = self.ms.total2matric(rx_)
        return np.append(collar_pot, rx)

    def solve_neumann(self, sim_time:float, t_act:list, rsx, cells:bool):
        """ solves the flux equations, with a neumann boundary condtion, see solve()
            @param sim_time [day]       needed for age dependent conductivities (age = sim_time - segment creation time)
            @param t_act [cm3 day-1]    tranpirational flux is negative
            @param rsx [cm]             soil matric potentials given per segment or per soil cell
            @param cells                indicates if the matric potentials are given per cell (True) or by segments (False)              
            @return [cm] root matric potential per root system node      
        """
        self.last = "neumann"
        if cells:
            rsx = self.get_hs(rsx)  # matric potential per root segment
        if not self.usecached_:
            self.update(sim_time)
        collar_pot = self.get_collar_potential(t_act, rsx)
        # print("solve_neumann(), collar potential", collar_pot, "cm", self.collar_index_)
        rsx_ = self.ms.matric2total(rsx)
        b = self.Kr.dot(rsx_)
        b[self.collar_index_] += self.kx0 * collar_pot
        rx = self.ms.total2matric(self.A_d_splu.solve(b))
        return np.append(collar_pot, rx)

    def solve(self, sim_time:float, t_act:list, rsx, cells:bool):
        """ Solves the hydraulic model using Neumann boundary conditions and switching to Dirichlet in case wilting point is reached
            @param sim_time [day]        needed for age dependent conductivities (age = sim_time - segment creation time)
            @param t_act [cm3 day-1]     transpiration rate
            @param rsx [cm]              soil matric potentials given per segment or per soil cell
            @param cells                 indicates if the matric potentials are given per cell (True) or by segments (False)
            @return [cm] root matric potential per root system node  
        """
        self.update(sim_time)
        #print('suf', self.suf.shape)
        if cells:
            rsx = self.get_hs(rsx)  # matric potential per root segment
        collar = self.get_collar_potential(t_act, rsx)
        collar = max(collar, self.wilting_point)
        rsx_ = self.ms.matric2total(rsx)
        b = self.Kr.dot(rsx_)
        b[self.collar_index_] += self.kx0 * collar
        rx = self.ms.total2matric(self.A_d_splu.solve(b))
        return np.append(collar, rx)

    def solve_again(self, sim_time:float, t_act:list, rsx, cells:bool):
        """ Solves the hydraulic model using Neumann boundary conditions and switching to Dirichlet in case wilting point is reached
            Depending of the method solve_again() is much faster using chached factorization of the last solve() command
            @param sim_time [day]        needed for age dependent conductivities (age = sim_time - segment creation time)
            @param t_act [cm3 day-1]     transpiration rate
            @param rsx [cm]              soil matric potentials given per segment or per soil cell
            @param cells                 indicates if the matric potentials are given per cell (True) or by segments (False)
            @return [cm] root matric potential per root system node  
        """
        if not self.cached:
            raise "HydraulicModel_Doussan.solve_again() makes only sense if cached = True"
            # self.update(sim_time)
        if cells:
            rsx = self.get_hs(rsx)  # matric potential per root segment
        collar = self.get_collar_potential(t_act, rsx)
        collar = max(collar, self.wilting_point)
        rsx_ = self.ms.matric2total(rsx)
        b = self.Kr.dot(rsx_)
        b[self.collar_index_] += self.kx0 * collar
        rx = self.ms.total2matric(self.A_d_splu.solve(b))
        return np.append(collar, rx)

    # def get_transpiration(self, sim_time, rx, rsx, cells = False):
    #     """ actual transpiration [cm3 day-1], calculated as the sum of all radial fluxes"""
    #     return np.sum(self.radial_fluxes(sim_time, rx, rsx, cells))

    def radial_fluxes(self, sim_time, rx, rsx, cells = False):
        """ returns the radial fluxes [cm3 day-1]"""
        if cells:
            rsx = self.get_hs(rsx)  # matric potential per root segment

        return -self.Kr.dot(rsx - rx[1:])  #   equals -q_root of Eqn (6) Leitner et al. (tba)

    def axial_fluxes(self, sim_time, rx, rsx = None, cells = None):
        """ returns the axial fluxes (independent of rsx in case of Doussan)
        @see axial_flux
        """
        n = self.ms.getNumberOfMappedSegments()
        return np.array([self.axial_flux(i, sim_time, rx) for i in range(0, n)])

    def axial_flux(self, seg_ind, sim_time, rx):
        """ returns the exact axial flux of segment ij of xylem model solution @param rx
            @param seg_ind              segment index
            @param sim_time [day]       needed for age dependent conductivities (age = sim_time - segment creation time)
            @param rx [cm]              root xylem matric potentials per root system node
            @return [cm3 day-1] axial volumetric flow rate
        """
        s = self.ms.segments[seg_ind]
        i, j = s.x, s.y
        n1, n2 = self.ms.nodes[i], self.ms.nodes[j]  # nodes
        v = n2.minus(n1)
        l = v.length()
        a = self.ms.getEffectiveRadius(seg_ind)  # radius
        st = int(self.ms.subTypes[seg_ind])  # sub type
        age = sim_time - self.ms.nodeCTs[int(s.y)]
        kr = self.params.kr_f(0, age, st, self.ot_root)  # c++ conductivity call back functions
        kx = self.params.kx_f(0, age, st, self.ot_root)  # c++ conductivity call back functi
        dpdz0 = (rx[j] - rx[i]) / l
        f = -kx * (dpdz0 - 1)
        return f

# int si, double age, int type, int orgtype

    def doussan_system_matrix(self, sim_time):
        """ """
        # print("doussan_system_matrix")
        IM = MappedPlantPython(self.ms).get_incidence_matrix()
        IMt = IM.transpose()
        kx_ = np.divide(self.params.getKx(sim_time), self.ms.segLength())  # / dl
        Kx = sparse.diags(kx_)
        kr = np.array(self.params.getEffKr(sim_time))
        # kr = np.maximum(np.ones(kr.shape) * 1.e-12, kr)
        Kr = sparse.diags(kr)
        L = IMt @ Kx @ IM  # Laplacian
        L_ = L[1:, 1:].tocsc()
        return  L_ + Kr, Kr, kx_[self.collar_index_]  # L_{N-1} + Kr, se Hess paper

    def update(self, sim_time):
        """ call before solve(), get_collar_potential(), and get_Heff() """
        # self.collar_index_ = self.collar_index()  # segment index of the collar segment
        A_d, self.Kr, self.kx0 = self.doussan_system_matrix(sim_time)
        self.A_d_splu = LA.splu(A_d)
        self.krs, _ = self.get_krs_(sim_time)
        self.suf = np.transpose(self.get_suf_())

    def get_collar_potential(self, t_act, rsx):
        """ collar potential for an actual transpiration (call update() before) """
        return (self.krs * self.get_heff_(rsx) - (-t_act)) / self.krs

    def get_soil_rootsystem_conductance(self, sim_time, h_bs, h_sr, sp):  # Vanderborgth et al. (2023), Eqn (12)
        """ The soil root system conductance per soil layer [day-1]
        
        sim_time             simulation time in days
        h_bs           bulk soil matric potential
        h_sr           matric potential at the soil root interface   
        sp             soil parameter: van Genuchten parameter set (type vg.Parameters)           
        """
        # krs, _ = self.get_krs(sim_time)  # [cm2/day]
        krs = self.krs
        # print("area", area)  #
        krs = krs / area  # [day-1]
        peri = Perirhizal(self.ms)  # helper class, wrap mappedSegments
        k_prhiz = peri.perirhizal_conductance_per_layer(h_bs, h_sr, sp)  # [day-1], Vanderborght et al. 2023, Eqn (6)
        # print("k_prhiz", np.nanmin(k_prhiz), np.nanmax(k_prhiz))
        # suf_ = self.get_suf(sim_time)  # [1]
        suf_ = self.suf
        suf = peri.aggregate(suf_[0,:])  # [1]
        # print("suf", np.min(suf), np.max(suf), np.sum(suf))
        return np.divide(krs * k_prhiz, suf * krs + k_prhiz)  # [day-1], see Vanderborgth et al. (2023), Eqn (12)

    def get_krs_(self, sim_time):
        """ calculatets root system conductivity [cm2/day] at simulation time @param sim_time [day] """
        # print("krs", sim_time)
        n = self.ms.getNumberOfMappedSegments()
        s = self.ms.segments[self.collar_index_]
        n2 = self.ms.nodes[s.y]
        rsx_ = np.ones((n, 1)) * (-500)  # total matric potential
        b = self.Kr.dot(rsx_)
        b[self.collar_index_, 0] += self.kx0 * -15000
        rx = self.A_d_splu.solve(b)  # total matric potential
        t_act = np.sum(-self.Kr.dot(rsx_ - rx))
        # print("get_krs() n2z", n2.z)
        # print("get_krs() rx[0]", rx[self.collar_index_, 0])
        # krs = -t_act / ((-500) - (rx[self.collar_index_, 0] - n2.z))
        krs = -t_act / ((-500) - (-15000))
        return krs, t_act

    def get_suf_(self):
        """ Standard uptake fraction (SUF) [1] per root segment, should add up to 1 """
        # print("suf")
        n = self.ms.getNumberOfMappedSegments()
        rsx = np.ones((n, 1)) * (-500)  # total matric potential
        b = self.Kr.dot(rsx)
        b[self.collar_index_, 0] += self.kx0 * -15000
        rx = self.A_d_splu.solve(b)  # total matric potential
        q = -self.Kr.dot(rsx - rx)
        return np.array(q) / np.sum(q)

    def get_heff_(self, rsx):
        """ effective total potential [cm] using cached suf """
        heff = self.suf.dot(self.ms.matric2total(rsx))
        # print("get_heff_()", heff[0], heff.shape)
        return heff[0]

