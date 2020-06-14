import timeit
import math

import numpy as np
from scipy import sparse
import scipy.sparse.linalg as LA

import plantbox as pb
from plantbox import XylemFlux
import rsml_reader as rsml


class XylemFluxPython(XylemFlux):
    """  Hybrid flux solver (following Meunier et al)"""

    def __init__(self, rs):
        """ @param rs is either a pb.MappedRootSystem or a string of a .rsml filename"""
        if isinstance(rs, str):
            rs = self.read_rsml(rs)
            super().__init__(rs)
        else:
            super().__init__(rs)

    def solve_neumann(self, sim_time :float, value :float, sxx) :
        """ solves the flux equations, with a neumann boundary condtion,
            @param sim_time [day]      simulation time to evaluate age dependent conductivities
            @param value [cm3 day-1]   tranpirational flux is negative 
         """
        # start = timeit.default_timer()
#         I, J, V, b = self.linear_system(sim_time)  # Python (care no age or type dependencies!)
#         self.aB = b
#         Q = sparse.coo_matrix((V, (I, J)))
        self.linearSystem(sim_time, sxx)  # C++
        Q = sparse.coo_matrix((np.array(self.aV), (np.array(self.aI), np.array(self.aJ))))
        Q = sparse.csr_matrix(Q)
        Q, b = self.bc_neumann(Q, self.aB, [0], [value / (self.rho * self.g)])  # cm3 day-1 -> something crazy (?)
        x = LA.spsolve(Q, b, use_umfpack = True)  # direct
        # print ("linear system assembled and solved in", timeit.default_timer() - start, " s")
        return x

    def solve_dirichlet(self, sim_time :float, value :float, sxc :float, sxx):
        """ solves the flux equations, with a dirichlet boundary condtion,
            @param sim_time [day]     simulation time to evaluate age dependent conductivities
            @param value [cm]         root collar pressure head 
            @param sx [cm]            soil pressure head around root collar segment 
         """
        self.linearSystem(sim_time, sxx)  # C++
        Q = sparse.coo_matrix((np.array(self.aV), (np.array(self.aI), np.array(self.aJ))))
        Q = sparse.csr_matrix(Q)
        Q, b = self.bc_dirichlet(Q, self.aB, [0], [float(value)])
        x = LA.spsolve(Q, b, use_umfpack = True)
        return x

    def solve(self, sim_time :float, value :float, sx :float, sxx, wilting_point :float = -15000):
        """ solves the flux equations using neumann and switching to dirichlet 
            in case wilting point is reached in root collar 
            @param simulation time  [day] for age dependent conductivities
            @param trans            [cm day-1] transpiration rate
            @param sx               [cm] soil matric potential around root collar
            @param sxx              [cm] soil matric potentials 
            @parm wiltingPoint      [cm] pressure head            
        """
        eps = 1.e-6
        x = [wilting_point - sx - 1]

        if sx >= wilting_point - eps:

            x = self.solve_neumann(sim_time, value, sxx)

            if x[0] < wilting_point - eps:
                Q = sparse.coo_matrix((np.array(self.aV), (np.array(self.aI), np.array(self.aJ))))
                Q = sparse.csr_matrix(Q)
                Q, b = self.bc_dirichlet(Q, self.aB, [0], [float(wilting_point)])
                try:
                    x = LA.spsolve(Q, b, use_umfpack = True)
                except:
                    print("Exeption solving Dirichlet")
                    print("Dirichlet at ", value - sx, "cm")
                    print("b", b)
                    raise
        else:
            print()
            print("solve_wp used Dirichlet because soil matric potential is below wilting point", sx)
            print()
            self.solve_dirichlet(sim_time, wilting_point, sx, sxx)

        return x

    def collar_flux(self, sim_time, rx, sx, seg_ind = 0):
        """ returns the exact transpirational flux of the solution @param rx [g/cm] """
        s = self.rs.segments[seg_ind]  # collar segment
        i, j = int(s.x), int(s.y)  # node indices
        n1, n2 = self.rs.nodes[i], self.rs.nodes[j]  # nodes
        v = n2.minus(n1)
        l = v.length()  # length of segment
        v.normalize()  # normalized v.z is needed for qz
        cell_ind = self.rs.seg2cell[seg_ind]
        p_s = sx[cell_ind]  # soil pressure at collar segment
        a = self.rs.radii[seg_ind]  # radius
        type = int(self.rs.types[seg_ind])  # conductivities kr, kx
        age = sim_time - int(self.rs.nodeCTs[j])
        kr = self.kr_f(age, type)  # c++ conductivity call back functions
        kx = self.kx_f(age, type)
        c = 2 * a * math.pi * kr / kx  # cm-2
        AA = np.array([[1, 1], [math.exp(math.sqrt(c) * l), math.exp(-math.sqrt(c) * l)] ])  # insert z = 0, z = l into exact solution
        bb = np.array([rx[i] - p_s, rx[j] - p_s])  # solve for solution
        d = np.linalg.solve(AA, bb)  # compute constants d_1 and d_2 from bc
        # p_r = lambda z: p_s + d[0] * math.exp(math.sqrt(c) * z) + d[1] * math.exp(-math.sqrt(c) * z)  # exact solution
        # dpdz = d[0] *sqrt(c)*exp(sqrt(c)*(-L)) + d[1] * (-sqrt(c)*exp(-sqrt(c)*(-L))) # dp/dz
        dpdz0 = d[0] * math.sqrt(c) - d[1] * math.sqrt(c)
        return -kx * (dpdz0 + v.z) * self.rho * self.g  # kx [cm5 day g-1]-> [cm3 day-1] by multiplying rho*g

    def get_nodes(self):
        """ converts the list of Vector3d to a 2D numpy array """
        return np.array(list(map(lambda x: np.array(x), self.rs.nodes)))

    @staticmethod
    def read_rsml(file_name :str):
        """ Reads and converts rsml to expected types and units"""
        polylines, props, funcs = rsml.read_rsml(file_name)
        nodes, segs = rsml.get_segments(polylines, props)
        radii, seg_ct, types = rsml.get_parameter(polylines, funcs, props)
        print("Read rsml:", len(nodes), "nodes", len(radii), "radii")
        nodes = np.array(nodes)  # for slicing in the plots
        nodes2 = []  # Conversions...
        for n in nodes:
            nodes2.append(pb.Vector3d(n[0] / 10., n[1] / 10., n[2] / 10.))  # [mm] -> [cm], and convert to a list of Vetor3d
        segs2 = []
        nodeCTs = np.zeros((len(nodes), 1))  # we need node creation times
        for i, s in enumerate(segs):
            nodeCTs[s[1]] = seg_ct[i]
            segs2.append(pb.Vector2i(int(s[0]), int(s[1])))
        radii = np.array(radii) / 10.  # [mm]->[cm]
        types = np.array(types, dtype = np.int64) - 1  # index must start with 0
        return pb.MappedSegments(nodes2, nodeCTs, segs2, radii, types)  # root system grid

    @staticmethod
    def bc_dirichlet(Q, b, n0, d):
        c = 0
        for c in range(0, len(n0)):
            i = n0[c]
            e0 = np.zeros((1, Q.shape[1]))  # build zero vector
            Q[i, :] = sparse.csr_matrix(e0)  # replace row i with ei
            Q[i, i] = 1
            b[i] = d[c]
        return Q, b

    @staticmethod
    def bc_neumann(Q, b, n0, f):
        c = 0
        for c in range(0, len(n0)):
            i = n0[c]  # print("Neumann BC at node "+str(i))
            b[i] += f[c]
        return Q, b

    @staticmethod
    def convert_(x, dtype = np.float64):
        return np.array(list(map(lambda x: np.array(x, dtype), x)), dtype)  # is there a better way?

    def linear_system(self, simTime :float):
        """ assembles the linear system (for comparison with the c++ equivalent linearSystem)"""
        Ns = len(self.rs.segments)
        N = len(self.rs.nodes)

        I = np.zeros(4 * Ns, dtype = np.int64)
        J = np.zeros(4 * Ns, dtype = np.int64)
        V = np.zeros(4 * Ns)
        b = np.zeros(N)
        k = 0
        for s in range(0, Ns):

            i, j = self.rs.segments[s].x, self.rs.segments[s].y

            a = self.rs.radii[j - 1]

            kx = self.kx[0]
            kr = self.kr[0]

            n1 = np.array([self.rs.nodes[i].x, self.rs.nodes[i].y, self.rs.nodes[i].z])
            n2 = np.array([self.rs.nodes[j].x, self.rs.nodes[j].y, self.rs.nodes[j].z])
            v = n2 - n1
            l = np.linalg.norm(v)
            vz = v[2] / l  # normed direction

            c = 2.*a * math.pi * kr / kx  # Eqn (2)
            d = math.exp(-math.sqrt(c) * l) - math.exp(math.sqrt(c) * l)  # Eqn (5)
            # print(i, j, n1, n2, l, d)
            di = 1. / d

            cii = -kx * di * math.sqrt(c) * (math.exp(-math.sqrt(c) * l) + math.exp(math.sqrt(c) * l))  # Eqn 16
            cij = 2 * kx * di * math.sqrt(c)  # Eqn 17
            bi = kx * vz  # Eqn 18 * self.rho * self.g

            b[i] += bi
            I[k], J[k], V[k] = i, i, cii
            k += 1
            I[k], J[k], V[k] = i, j, cij
            k += 1

            i, j = j, i  # edge ji
            b[i] -= bi  # Eqn 14 with changed sign
            I[k], J[k], V[k] = i, i, cii
            k += 1
            I[k], J[k], V[k] = i, j, cij
            k += 1

        return I, J, V, b

