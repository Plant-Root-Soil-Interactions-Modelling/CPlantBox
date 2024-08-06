import sys; sys.path.append(".."); sys.path.append("../..");

import timeit

import numpy as np
from scipy import sparse
import scipy.sparse.linalg as LA
import matplotlib.pyplot as plt

import plantbox as pb
from plantbox import XylemFlux
import rsml.rsml_reader as rsml


def sinusoidal(t):
    """ sinusoidal function (used for transpiration) (integral over one day is 1)"""
    return np.sin(2. * np.pi * np.array(t) - 0.5 * np.pi) + 1.


def sinusoidal2(t, dt):
    """ sinusoidal function from 6:00 - 18:00, 0 otherwise (integral over one day is 1)"""
    return np.maximum(0., np.pi * (np.cos(2 * np.pi * (t - 0.5)) + np.cos(2 * np.pi * ((t + dt) - 0.5))) / 2)


class XylemFluxPython(XylemFlux):
    """  Hybrid flux solver (following Meunier et al.)
    
        C++ part is defined in CPlantBox XylemFlux.hh and XylemFlux.cpp
        
        Calculates water movement within the xylems, assuming a constant matric potential around each xylem segment,
        given for each segment, or for each soil cell. 
        
        The root surface flux is calculated exactly as in Meunier et al.        
    """

    def __init__(self, rs):
        """ @param rs is either a pb.MappedRootSystem, pb.MappedSegments, or a string containing a rsml filename"""
        if isinstance(rs, str):
            rs = self.read_rsml(rs)
            super().__init__(rs)
        else:
            super().__init__(rs)

        self.neumann_ind = [0]  # node indices for Neumann flux
        self.dirichlet_ind = [0]  # node indices for Dirichlet flux
        self.collar_index_ = self.collar_index()  # segment index of the collar segement

        self.last = "none"
        self.Q = None  # store linear system
        self.b = None

        self.Vmax = 45.25 * 62 * 1.e-11 * (24.*3600.)  # kg/(m2 day) -> York et. al (2016) (Lynch group)
        self.Km = 10.67 * 62 * 1.e-6  # kg/m3
        self.CMin = 4.4 * 62 * 1.e-6  # kg/m3
        self.Exu = 1.e-2  # data Eva Oburger kg /(m2 day)
        self.plant = rs

    def solute_fluxes(self, c):
        """ concentrations @param c [kg/m3], returns [g/day]
        self.Vmax [kg/(m2 day)]
        self.Km [kg/m3]
        self.Cmin [kg/m3]
        """
        segs = self.rs.segments
        a = self.rs.radii
        l = self.rs.segLength()
        assert(len(segs) == len(c))
        sf = np.zeros(len(segs),)
        c = np.maximum(c, self.CMin)  ###############################################
        for i, s in enumerate(segs):
            sf[i] = -2 * np.pi * a[i] * l[i] * 1.e-4 * ((c[i] - self.CMin) * self.Vmax / (self.Km + (c[i] - self.CMin)))  # kg/day
        sf = np.minimum(sf, 0.)
        return sf * 1.e3  # kg/day -> g/day

    def exudate_fluxes(self, rs_age, kex):
        """ returns [g/day]
        self.Exu [kg/(m2 day)]
        """
        segs = self.rs.segments
        tipI = np.array(self.rs.getRootTips()) - 1
        polylengths = self.rs.getParameter("length")
        types = self.rs.getParameter("type")

        a = self.rs.radii
        l = self.rs.segLength()
        sf = np.zeros((len(segs)))
        kex_all = np.zeros((len(segs)))

        for i, s in enumerate(tipI):
            if types[i] != 0:
                l_ = 0
                s_ = s
                while l_ <= kex[0][1]:
                    l_ = l_ + l[s_]
                    kexu = (kex[1][1] - kex[1][0]) / (kex[0][1] - kex[0][0]) * l_ + kex[1][0]  # linear decrease
                    kex_all[s_] = max(0, kexu)
                    sf[s_] = 2 * np.pi * a[s_] * l[s_] * 1.e-4 * kex_all[s_] * 1.e3  # g/day/root
                    s_ = s_ - 1
                    if l_ > polylengths[i]:
                        break

        return sf, kex_all  #  g/day

    def get_incidence_matrix(self):
        """ retruns the incidence matrix (number of segments, number of nodes) of the root system in self.rs 
        """
        segs = self.rs.segments
        sn = len(segs)
        nn = len(self.rs.nodes)
        ii_, jj_, vv_ = [], [], []
        for i, s in enumerate(segs):  # build incidence matrix from edges
            ii_.append(i)
            ii_.append(i)
            jj_.append(segs[i].x)
            jj_.append(segs[i].y)
            vv_.append(-1.)
            vv_.append(1.)
        return sparse.coo_matrix((np.array(vv_), (np.array(ii_), np.array(jj_))), shape = (sn, nn))

    def linearSystem_doussan(self, sim_time, sxx, cells = True, soil_k = []):
        """ soil_k TODO
        """
        print("building doussan matrix")
        IM = self.get_incidence_matrix()
        IMt = IM.transpose()
        kx_ = np.divide(self.getKx(sim_time), self.rs.segLength())  # / dl
        Kx = sparse.diags(kx_)
        kr = np.zeros((IM.shape[1],))
        kr[1:] = np.array(self.getEffKr(sim_time))  #  change to getKr(), use that instead
        Kr = sparse.diags(kr)
        A = IMt @ Kx @ IM + Kr
        (self.aI, self.aJ, self.aV) = sparse.find(A)
        if cells:
            hs = np.zeros((IM.shape[1],))
            hs[1:] = self.getKs(sim_time)
            self.aB = Kr * hs
        else:
            hs = np.zeros((IM.shape[1],))
            hs[1:] = sxx
            self.aB = Kr * hs

    def solve_neumann(self, sim_time:float, value, sxx, cells:bool, soil_k = []):
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
            n = len(self.neumann_ind)
            value = [value / n] * n

        if len(soil_k) > 0:
            self.linearSystem(sim_time, sxx, cells, soil_k)  # C++ (see XylemFlux.cpp)
        else:
            self.linearSystem(sim_time, sxx, cells)  # C++ (see XylemFlux.cpp)

        self.Q = sparse.coo_matrix((np.array(self.aV), (np.array(self.aI), np.array(self.aJ))))
        self.Q = sparse.csc_matrix(self.Q)
        self.Q, self.b = self.bc_neumann(self.Q, self.aB, self.neumann_ind, value)  # cm3 day-1

        x = LA.spsolve(self.Q, self.b, use_umfpack = True)  # direct
        return x

    def solve_dirichlet(self, sim_time:float, value:list, sxc:float, sxx, cells:bool, soil_k = []):
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
            n = len(self.dirichlet_ind)
            value = [value] * n

        if len(soil_k) > 0:
            self.linearSystem(sim_time, sxx, cells, soil_k)  # C++ (see XylemFlux.cpp)
        else:
            self.linearSystem(sim_time, sxx, cells)

        self.Q = sparse.coo_matrix((np.array(self.aV), (np.array(self.aI), np.array(self.aJ))))
        self.Q = sparse.csc_matrix(self.Q)
        self.Q, self.b = self.bc_dirichlet(self.Q, self.aB, self.dirichlet_ind, value)

        x = LA.spsolve(self.Q, self.b, use_umfpack = True)
        return x

    def solve(self, sim_time:float, trans:list, sx:float, sxx, cells:bool, wilting_point:float, soil_k = []):
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
            self.last = "neumann"

            if x[0] <= wilting_point:

                Q = sparse.coo_matrix((np.array(self.aV), (np.array(self.aI), np.array(self.aJ))))
                Q = sparse.csc_matrix(Q)
                Q, b = self.bc_dirichlet(Q, self.aB, [0], [float(wilting_point)])
                x = LA.spsolve(Q, b, use_umfpack = True)
                self.last = "dirichlet"

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
        s = self.rs.segments[seg_ind]
        nodes = self.rs.nodes
        numleaf = 0
        organTypes = self.get_organ_types()
        ot = int(organTypes[seg_ind])  # for conductivities kr, kx
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
        kx = self.kx_f(age, st, ot, seg_ind)
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

    def collar_flux(self, sim_time, rx, sxx, k_soil = [], cells = True):
        """ returns the exact transpirational flux of the xylem model solution @param rx
            @see axial_flux        
        """
        self.collar_index_ = self.collar_index()
        return self.axial_flux(self.collar_index_, sim_time, rx, sxx, k_soil, cells, ij = True)

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
        return np.array(self.segFluxes(sim_time, rx, sxx, False, cells, k_soil))  # approx = False

    def get_nodes(self):
        """ converts the list of Vector3d to a 2D numpy array """
        return np.array(list(map(lambda x: np.array(x), self.rs.nodes)))

    def get_segments(self):
        """ converts the list of Vector2i to a 2D numpy array """
        return np.array(list(map(lambda x: np.array(x), self.rs.segments)), dtype = np.int64)

    def get_subtypes(self):
        """ segment sub types as numpy array """
        return np.array(self.rs.subTypes)

    def get_ages(self, final_age = -1.):
        """ converts the list of nodeCT to a numpy array of segment ages
        @param final_age [day]         current root system age, (default = 0 means detect maximum from nodeCT)
        """
        cts = np.array(self.rs.nodeCTs)
        if final_age == -1.:
            final_age = np.max(cts)
        node_ages = final_age * np.ones(cts.shape) - cts  # from creation time to age
        segs = self.get_segments()
        ages = np.zeros(segs.shape[0])
        for i, s in enumerate(segs):
            ages[i] = node_ages[s[1]]  # segment age based on second node (as in XylemFlux::linearSystem)
        return ages

    def get_conductivities(self, final_age = -1.):
        """ returns radial [1 day-1] and axial [cm3 day-1] conductivity per segment 
        @param final_age [day]         current root system age, (default = 0 means detect maximum from nodeCT)        
        """
        ages = self.get_ages(final_age)
        subtypes = self.get_subtypes()
        assert ages.shape == subtypes.shape, "get_conductivities: length of age and subtype must agree"
        kr = np.zeros(ages.shape)
        kx = np.zeros(ages.shape)
        for i in range(0, len(ages)):
            kr[i] = self.kr_f(ages[i], subtypes[i])
            kx[i] = self.kx_f(ages[i], subtypes[i])
        return kr, kx

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
        """ return segment indices of organ type @param ot """
        organTypes = self.get_organ_types()
        segIdx = np.array(list(range(0, len(organTypes))))
        otsegs = segIdx[organTypes == ot]
        return otsegs

    def get_organ_types(self):
        """ segment organ types as numpy array """
        return np.array(self.rs.organTypes)

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

    def get_suf(self, sim_time, approx = False, organType_ = None):
        """ calculates the standard uptake fraction [1] at simulation time @param sim_time [day]
            (suf is constant for age independent conductivities) 
            for the root system  (organType_=2 or organType_=None)
            the shoot (organType_=4, @see get_krs(plant = True)) 
        """
        segs = self.rs.segments
        nodes = self.rs.nodes
        p_s = np.zeros((len(segs),))
        organType = self.get_organ_types()
        for i, s in enumerate(segs):
            if (organType[i] == 2):
                p_s[i] = -500 - 0.5 * (nodes[s.x].z + nodes[s.y].z)  # constant total potential (hydraulic equilibrium)
            else:
                p_s[i] = self.airPressure - 0.5 * (nodes[s.x].z + nodes[s.y].z)  # constant total potential (hydraulic equilibrium)
        if(pb.leaf in organType):  # we have a whole plant
            rx = self.solve_neumann(sim_time, 0, p_s, cells = False)
            # transpiration == sum(leaf radial flux) == sum(root radial flux)
        else:
            transpiration = -1.e5
            rx = self.solve_neumann(sim_time, transpiration , p_s, cells = False)  # False: matric potential not given per cell (but per segment)

        fluxes = np.array(self.segFluxes(sim_time, rx, p_s, approx = approx, cells = False))  # cm3/day, simTime,  rx,  sx,  approx, cells

        if(pb.leaf in organType):
            if organType_ == pb.leaf:  # we want the suf for the leaves, only usefull when computing get_krs() for plants
                transpiration = sum(np.array(fluxes)[organType == pb.leaf])
            else:  # we want the suf for the roots
                transpiration = -sum(np.array(fluxes)[organType == pb.leaf])

        if not (organType_ is None):  # we want the suf for only specific orgran types
            fluxes = fluxes[organType == organType_]
        else:
            fluxes[organType == pb.leaf] = np.nan

        if transpiration != 0:
            return fluxes / transpiration  # [1]

        else:
            return np.full(len(np.array(fluxes)[organType == organType_]), 0.)

    def get_mean_suf_depth(self, sim_time):
        """  mean depth [cm] of water uptake based suf """
        suf = self.get_suf(sim_time)
        suf[np.isnan(suf)] = 0  # to get the mean in spite of the nan
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

    def get_krs(self, sim_time, seg_ind = [0], plant = False):
        """ calculatets root system conductivity [cm2/day] (plant=False)
        or root, shoot and total conductivity (plant=True)
        at simulation time @param sim_time [day] 
        if there is no single collar segment at index 0, pass indices using @param seg_ind, see find_base_segments        
        """
        segs = self.rs.segments
        nodes = self.rs.nodes
        p_s = np.zeros((len(segs),))
        organType = self.get_organ_types()
        subType = self.get_subtypes()
        for i, s in enumerate(segs):
            if (organType[i] == 2):
                p_s[i] = -500 - 0.5 * (nodes[s.x].z + nodes[s.y].z)  # constant total potential (hydraulic equilibrium)
            else:
                p_s[i] = self.airPressure - 0.5 * (nodes[s.x].z + nodes[s.y].z)  # constant total potential (hydraulic equilibrium)

        if plant:
            rx = self.solve_neumann(sim_time, 0., p_s, cells = False)
            fluxes = self.segFluxes(sim_time, rx, p_s, approx = False, cells = False)  # cm3/day, simTime,  rx,  sx,  approx, cells
            jc = sum(np.array(fluxes)[organType == 4])
            root_p_s = p_s[organType == 2]
            eswp = self.get_eswp(sim_time, root_p_s, cells = False, organType_ = 2)
            leaf_p_s = p_s[organType == 4]

            eawp = self.get_eswp(sim_time, leaf_p_s, cells = False, organType_ = 4)

            rootOs = self.rs.getOrgans(2)
            ShootRoots = [rootO for rootO in rootOs if rootO.getParameter("subType") == self.rs.getSeed().shootborneType]  # [0]
            if len(ShootRoots) > 0:  # ATT! does not work if tiller separation is below the higher Shoot-root
                ShootRoot = ShootRoots[0]
                rootCollar = ShootRoot.getNodeId(0)
            else:  # noo shoot root: we can take the wat pot at the seed
                rootCollar = 0

            if eswp == 0:  # no transpiration
                print("esp == 0, no transpiration?")
                return np.nan, np.nan, np.nan, 0, np.nan, np.nan, rx[rootCollar]

            kr_roots = jc / (eswp - rx[rootCollar])
            kr_shoot = jc / (rx[rootCollar] - eawp)
            kr_plant = jc / (eswp - eawp)
            # conductivites, transpirations and equvalient water potentials
            return kr_roots, kr_shoot, kr_plant , jc, eswp, eawp, rx[rootCollar]
        else:
            rx = self.solve_dirichlet(sim_time, -15000, 0., p_s, cells = False)  # p_s is a constant total potential opf -500 cm
            jc = 0
            for i in seg_ind:
                jc -= self.axial_flux(i, sim_time, rx, p_s, [], cells = False, ij = True)  # Transpiration jc according to Meunier et al.
            # jc = sum Q = sum Krs SUF (Heff - Hcollar) = Krs (Heff - Hcollar)
            # krs = jc / (Heff - Hcollar) # see Eqn (8) in Vanderborght et al (2021), it is more general then using C4 as in Eqn (14)
            krs = jc / (-500 - 0.5 * (nodes[segs[0].x].z + nodes[segs[0].y].z) - rx[self.dirichlet_ind[0]])
            # Hcollar =  rx[self.dirichlet_ind[0] + 0.5 * (nodes[segs[0].x].z + nodes[segs[0].y].z) # matric -> total potential
            return krs , jc

    def get_eswp(self, sim_time, p_s, cells = True, organType_ = None):
        """ calculates the equivalent soil water potential [cm] (organType_ = 2 or organType_ = None)
        or air water potential (organType_ = 4) (@see get_krs(plant = True))
        at simulation time @param sim_time [day] for 
        the potential @param p_s [cm] given per cell """
        organType = self.get_organ_types()

        segs = self.get_segments()
        nodes = self.rs.nodes
        seg2cell = self.rs.seg2cell
        suf = self.get_suf(sim_time, organType_ = organType_)
        if not (organType_ is None):
            segs = segs[organType == organType_]  # only happens/useful when using get_krs() for plants
        else:
            suf[np.isnan(suf)] = 0  # to get the mean in spite of the nan
        eswp = 0.
        for i, s in enumerate(segs):
            if cells:
                eswp += suf[i] * (p_s[seg2cell[i]] + 0.5 * (nodes[s[0]].z + nodes[s[1]].z))  # matric potential to total potential
            else:
                eswp += suf[i] * (p_s[i] + 0.5 * (nodes[s[0]].z + nodes[s[1]].z))  # matric potential to total potential
        return eswp

    def kr_f(self, age, st, ot = 2 , seg_ind = 0, cells = False):
        """ root radial conductivity [1 day-1] for backwards compatibility """
        return self.kr_f_cpp(seg_ind, age, st, ot, cells)  # kr_f_cpp is XylemFlux::kr_f

    def kx_f(self, age, st, ot = 2, seg_ind = 0):
        """ root axial conductivity [cm3 day-1]  for backwards compatibility """
        return self.kx_f_cpp(seg_ind, age, st, ot)  # kx_f_cpp is XylemFlux::kx_f

    def collar_index(self):
        """ returns the segment index of the collar segment """
        segs = self.rs.segments
        for i, s in enumerate(segs):
            if s.x == 0:
                return i

    def test(self):
        """ perfoms some sanity checks, and prints to the console """
        print("\nXylemFluxPython.test():")
        # 1 check if segment index is node index-1
        segments = self.get_segments()
        nodes = self.get_nodes()
        types = self.rs.subTypes
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
        seg_length = self.rs.segLength()
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
        map = self.rs.seg2cell
        for seg_id, cell_id in map.items():
            if cell_id < 0:
                print("Warning: segment ", seg_id, "is not mapped, this will cause problems with coupling!", nodes[segments[seg_id][0]], nodes[segments[seg_id][1]])
        print()

    def plot_conductivities(self, monocot = True, plot_now = True, axes_ind = [], lateral_ind = []):
        """ plots conductivity  
        @param monocot      indicates if monocot (True) or dicot (False)
        @param plot_now     indicates if the figure is shown, or just retruned
        @param axes_ind     for monocots a list of three root types for "tap root", "basal", "shoot borne";
                            for dicots a list of one root type representing the "tap root"
        @param lateral_ind  for monocots a list of two root types for "1st order laterals", "2nd order laterals"
                            for dicots a list of three root types for "1st order laterals", "2nd order laterals", "3rd order laterals"          
        """
        axes_age = np.linspace(-5, 100, 500)
        lateral_age = np.linspace(-2, 25, 125)
        lateral_cols = ["r", "g:", "m--", "b--"]
        axes_cols = ["r", "g:", "m--", "b--"]
        if monocot:
            axes_str = ["tap root", "basal", "shoot borne"]
            lateral_str = ["1st order laterals", "2nd order laterals"]
            if axes_ind == []:
                axes_ind = [1, 4, 5]
            if lateral_ind == []:
                lateral_ind = [2, 3]
        else:  # dicot
            axes_str = ["tap root"]
            lateral_str = ["1st order laterals", "2nd order laterals", "3rd order laterals"]
            if axes_ind == []:
                axes_ind = [1]
            if lateral_ind == []:
                lateral_ind = [2, 3, 4]

        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize = (16, 10))
        for j, st in enumerate(axes_ind):
            kx_ = [ self.kx_f(axes_age[i], st, 2) for i in range(0, len(axes_age)) ]
            ax1.plot(axes_age, kx_, axes_cols[j])
        kx_max = np.max(kx_)
        ax1.legend(axes_str)
        ax1.set_title("Axis")
        ax1.set_xlabel("age [day]")
        ax1.set_ylabel("axial conductance [cm$^3$ day$^{-1}$]")
        for j, st in enumerate(lateral_ind):
            kx_ = [ self.kx_f(lateral_age[i], st, 2) for i in range(0, len(lateral_age)) ]
            ax2.plot(lateral_age, kx_, axes_cols[j])
        kx_max = max(kx_max, np.max(kx_))
        ax2.legend(lateral_str)
        ax2.set_title("Laterals")
        ax2.set_xlabel("age [day]")
        ax2.set_ylabel("axial conductance [cm$^3$ day$^{-1}$]")
        for j, st in enumerate(axes_ind):
            kr_ = [ self.kr_f(axes_age[i], st, 2, 0) for i in range(0, len(axes_age)) ]
            ax3.plot(axes_age, kr_, axes_cols[j])
        kr_max = np.max(kr_)
        ax3.legend(axes_str)
        ax3.set_title("Axis")
        ax3.set_xlabel("age [day]")
        ax3.set_ylabel("radial conductance [day$^{-1}$]")
        for j, st in enumerate(lateral_ind):
            kr_ = [ self.kr_f(lateral_age[i], st, 2, 0) for i in range(0, len(lateral_age)) ]
            ax4.plot(lateral_age, kr_, axes_cols[j])
        kr_max = max(kr_max, np.max(kr_))
        ax4.legend(lateral_str)
        ax4.set_title("Laterals")
        ax4.set_xlabel("age [day]")
        ax4.set_ylabel("radial conductance [day$^{-1}$]")
        print(kx_max)
        print(kr_max)
        ax1.set_ylim([0, kx_max * 1.1])
        ax2.set_ylim([0, kx_max * 1.1])
        ax3.set_ylim([0, kr_max * 1.1])
        ax4.set_ylim([0, kr_max * 1.1])
        print()
        for st in range(0, 5):
            print("SubType {:g} for negative age: kx = {:g}, kr = {:g}".format(st, self.kx_f(-1, st, 2), self.kr_f(-1, st, 2, 0)))
            print("SubType {:g} for old root age: kx = {:g}, kr = {:g}".format(st, self.kx_f(100, st, 2), self.kr_f(100, st, 2, 0)))
        print()
        if plot_now:
            plt.show()
        return fig

    def summarize_fluxes(self, fluxes, sim_time, rx, p_s, k_soil = [], cells = False, show_matrices = False):
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
        axialfluxes_j = np.array([self.axial_flux(i, sim_time, rx = rx, sxx = p_s, k_soil = k_soil, ij = True, cells = cells) for i in range(0, len(self.rs.segments))])
        axialfluxes_i = -np.array([self.axial_flux(i, sim_time, rx = rx, sxx = p_s, k_soil = k_soil, ij = False, cells = cells) for i in range(0, len(self.rs.segments))])

        balance = -axialfluxes_i - axialfluxes_j + fluxes
        if show_matrices:
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
        get_y_node = lambda vec: np.array(vec)[1]
        get_x_node = lambda vec: np.array(vec)[0]
        get_nodetype = lambda y: organTypes[y - 1]
        nodesy = np.array([get_y_node(xi) for xi in segments])
        nodesx = np.array([get_x_node(xi) for xi in segments])
        nodeflux = np.full(len(self.get_nodes()), np.nan)
        nodes = np.arange(len(self.get_nodes()))
        for nd in range(len(nodeflux)):
            fluxseg_roots_bellow = axialfluxes_j[np.logical_and((nodesx == nd), (organTypes == 2))]  # upper root flux where node is upper node
            fluxseg_roots_above = axialfluxes_i[np.logical_and((nodesy == nd), (organTypes == 2))]  # lower root flux where node is lower node
            fluxseg_stemsleaves_bellow = axialfluxes_j[np.logical_and((nodesy == nd), (organTypes > 2))]  # upper root flux where node is upper node
            fluxseg_stemsleaves_above = axialfluxes_i[np.logical_and((nodesx == nd), (organTypes > 2))]  # lower root flux where node is lower node
            nodeflux[nd] = sum((sum(fluxseg_roots_bellow), sum(fluxseg_roots_above), sum(fluxseg_stemsleaves_bellow), sum(fluxseg_stemsleaves_above)))
        tiproots, tipstems, tipleaves = self.get_organ_nodes_tips()
        tips = np.concatenate((tiproots, tipstems, tipleaves))
        nodeflux[tips] = np.nan
        print('sum (absolut value) of net water balance in each node, appart from node at the tip of organs {} [cm3/day] (should be 0)'.format(np.nansum(abs(nodeflux))))
        tot_balance = -sum(axialfluxes_i[tipsegroots]) - sum(axialfluxes_j[tipsegstems]) - sum(axialfluxes_j[tipsegleaves]) - sum(fluxes)
        print('net water flux in plant: {} [cm3/day] (should be 0)'.format(tot_balance))
        return fluxes

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
            print()

        return pb.MappedSegments(nodes2, cts, segs2, segRadii, segTypes)  # root system grid

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
