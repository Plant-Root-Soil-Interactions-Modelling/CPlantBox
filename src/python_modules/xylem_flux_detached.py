from xylem_flux import *


class XylemFluxDetached(XylemFluxPython):
    
    def __init__(self, rs):
        """ @param rs is either a pb.MappedRootSystem, pb.MappedSegments, or a string containing a rsml filename"""
        if isinstance(rs, str):
            rs = self.read_rsml(rs)
            super().__init__(rs)
        else:
            super().__init__(rs)

        self.seg_ind = [0]  # segment indices for Neuman flux
        self.node_ind = [0]  # node indices for Dirichlet flux

    def solve_neumann(self, sim_time:float, value, sxx, cells:bool, soil_k=[]):
        """ solves the flux equations, with a neumann boundary condtion, see solve()
            @param sim_time [day]       needed for age dependent conductivities (age = sim_time - segment creation time)
            @param value [cm3 day-1]    tranpirational flux is negative
            @param sxx [cm]             soil matric potentials given per segment or per soil cell
            @param cells                indicates if the matric potentials are given per cell (True) or by segments (False)
            @param soil_k [day-1]       optionally, soil conductivities can be prescribed per segment, 
                                        conductivity at the root surface will be limited by the value, i.e. kr = min(kr_root, k_soil)  
            @return [cm] root xylem pressure per root system node         
         """
        print("neumann detached")
         
        # start = timeit.default_timer()
        if isinstance(value, (float, int)):
            n = len(self.seg_ind)
            value = [value / n] * n

        if len(soil_k) > 0:
            self.linearSystem_detached(sim_time, sxx, cells, soil_k)  # C++ (see XylemFlux.cpp)
        else:
            self.linearSystem_detached(sim_time, sxx, cells)  # C++ (see XylemFlux.cpp)

        Q = sparse.coo_matrix((np.array(self.aV), (np.array(self.aI), np.array(self.aJ))))
        Q = sparse.csr_matrix(Q)
        Q, b = self.bc_neumann(Q, self.aB, self.seg_ind, value)  # cm3 day-1
        x = LA.spsolve(Q, b, use_umfpack=True)  # direct
        # print ("linear system assembled and solved in", timeit.default_timer() - start, " s")
        return x

    def solve_dirichlet(self, sim_time:float, value:list, sxc:float, sxx, cells:bool, soil_k=[]):
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
        print("dirichlet detached") 
         
        if isinstance(value, (float, int)):
            n = len(self.node_ind)
            value = [value] * n

        if len(soil_k) > 0:
            self.linearSystem_detached(sim_time, sxx, cells, soil_k)  # C++ (see XylemFlux.cpp)
        else:
            self.linearSystem_detached(sim_time, sxx, cells)

        Q = sparse.coo_matrix((np.array(self.aV), (np.array(self.aI), np.array(self.aJ))))
        Q = sparse.csr_matrix(Q)
        Q, b = self.bc_dirichlet(Q, self.aB, self.node_ind, value)
        x = LA.spsolve(Q, b, use_umfpack=True)
        return x
