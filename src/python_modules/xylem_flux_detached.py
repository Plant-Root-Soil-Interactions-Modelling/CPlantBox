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
#        print("neumann detached")

        # start = timeit.default_timer()
        if isinstance(value, (float, int)):
            n = len(self.seg_ind)
            value = [value / n] * n

        if len(soil_k) > 0:
            self.linearSystem_detached(sim_time, sxx, cells, soil_k)  # C++ (see XylemFlux.cpp)
        else:
            self.linearSystem_detached(sim_time, sxx, cells)  # C++ (see XylemFlux.cpp)

        Q = sparse.coo_matrix((np.array(self.aV), (np.array(self.aI), np.array(self.aJ))))
        Q = sparse.csc_matrix(Q)
        Q, b = self.bc_neumann(Q, self.aB, self.seg_ind, value)  # cm3 day-1
        x = LA.spsolve(Q, b, use_umfpack=True)  # direct

#         plt.spy(Q)
        Qinv = LA.inv(Q)
        # plt.spy(Qinv)
        # q = Qinv.toarray().flatten()
        # print(np.mean(q), np.std(q))
        # plt.hist(q.flatten(), bins = 100)

        plt.show()

        # plt.show()
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
        # print("dirichlet detached")

        if isinstance(value, (float, int)):
            n = len(self.node_ind)
            value = [value] * n

        if len(soil_k) > 0:
            self.linearSystem_detached(sim_time, sxx, cells, soil_k)  # C++ (see XylemFlux.cpp)
        else:
            self.linearSystem_detached(sim_time, sxx, cells)

        Q = sparse.coo_matrix((np.array(self.aV), (np.array(self.aI), np.array(self.aJ))))
        Q = sparse.csc_matrix(Q)
        Q, b = self.bc_dirichlet(Q, self.aB, self.node_ind, value)
        x = LA.spsolve(Q, b, use_umfpack=True)
        return x
    
    def mean_xylem_pressure(self, seg_ind, sim_time, rx, sxx, cells=True):
        s = self.rs.segments[seg_ind]        
        nodes = self.rs.nodes
        organTypes = self.get_organ_types()
        ot = int(organTypes[seg_ind])  # for conductivities kr, kx
        i = 0  # DETACHED !!!
        j = int(s.y)
        n1, n2 = self.rs.nodes[s.x], self.rs.nodes[j]  # nodes
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
        kr = self.kr_f(age, st, ot, 0, seg_ind)  # c++ conductivity call back functions
        kx = self.kx_f(age, st, ot, seg_ind)
        if a * kr > 1.e-16:
            tau = math.sqrt(2 * a * np.pi * kr / kx)  # cm-2
            AA = np.array([[1, 1], [np.exp(tau * l), np.exp(-tau * l)] ])
            bb = np.array([rx[i] - p_s, rx[j] - p_s])  # solve for solution
            d = np.linalg.solve(AA, bb)  # compute constants d_1 and d_2 from bc
            mean_p = (1. / l) * ((np.exp(-tau * l) * (np.exp(tau * l) - 1) * (d[0] * np.exp(tau * l) + d[1])) / tau + l * p_s)
        else:  # solution for kr = 0, or a = 0
            mean_p = (rx[j] + rx[i]) * 0.5

#         if j == 1:
#             print(rx[i], rx[j], p_s, mean_p, l, kr, kx, a)

        return mean_p
    
