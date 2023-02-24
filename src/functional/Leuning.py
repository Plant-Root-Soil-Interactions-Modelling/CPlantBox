
import numpy as np
from scipy import sparse
import scipy.sparse.linalg as LA
from functional.xylem_flux import XylemFluxPython  # Python hybrid solver
# from plantbox import PhloemFlux


class Leuning(XylemFluxPython):  # PhloemFlux,
    """      
        C++ part is defined in CPlantBox XylemFlux.hh and XylemFlux.cpp
        Calculates water movement within the xylems, the stomatal aperture and net assimilation rate. 
        
        The root surface flux is calculated exactly as in Meunier et al.        
    """

    def __init__(self, rs):
        """ @param rs is either a pb.MappedPlant or a string containing a rsml filename
            sets the parameters for the calculation of gs and of the transpiration rate
            from Tuzet 2003, Dewar 2002 and Lu 2020                     
        """
        super().__init__(rs)
        self.Param = {  # list of parameters in alphabetical order. TODO: delete the useless parameters
        'a1':3,  # 6, #(-) fitting parameter for gco2
        'a2': 1.6,  # (-)
        'a3':1.7,
        'alpha': 0.2,
        'alpha2': 1.5,
        'anc': 5.35,
        'arn': 1.8351,
        'beta': 5,
        'bnc': 0.44,
        'brn':-1.5075,
        'c': 7.72 * 5,  # added *5 with regards to Lu2020
        'Chl/N': 3.7,
        'cs':350e-6,  # example from Dewar2002
        'd': 0.04,
        'D0': 1.67,  # kPa
        'DMcrit': 1.55,
        'Eac': 59430,
        'Eaj': 37000,
        'Eao': 36000,
        'Eard': 53000,
        'Eav': 58520,
        'Edj': 220000,
        'Edv': 220000,
        'FRN': 0.28,
        'fwr': 9.308 * 1e-2,  # parametrised with data from corso2020
        'g0': 0.0,
        'gamma0': 28e-6,
        'gamma1': 0.0509,
        'gamma2': 0.001,
        'JmaxrefChl1':2.78,
        'JmaxrefChl2':18.45,
        'Kc_ref': 302e-6,
        'Ko_ref': 256e-3,
        'Mh2o': 18,
        'Ncmin': 4.4,
        'oi': 210e-3,
        'Patm': 101,  # kPa
        'pi_g0':-1000,  # kpa
        'p_lcrit':-0.869,
        'R': 8.3143,
        'Rd_ref': 0.32e-6,
        'rho_h2o': 1,
        'rho_p': 0.35,
        'S': 700,
        'sh': 3.765 * 1e-4,
        'Tref': 293.2,
        'theta': 0.9,
        'VcmaxrefChl1':1.28 / 2,  # otherwise value too high
        'VcmaxrefChl2':8.33 / 2,  # otherwise value too high
        }

    def solve_leuning(self, sim_time:float, sxx, cells:bool, Qlight , VPD: float, Tl, p_linit, ci_init, cs, soil_k = [], N = [], log = True, verbose = False):
        """ solves the flux equations, with a neumann boundary condtion, see solve()
            @param sim_time [day]           needed for age dependent conductivities (age = sim_time - segment creation time)
            @param sxx [cm]                 soil matric potentials given per segment or per soil cell
            @param cells                    indicates if the matric potentials are given per cell (True) or by segments (False)
            @param Q [mol photons m-2 s-1]  absorbed photon irradiance (mean or per leaf segment)
            @param Tl [K]                   leaf temperature (mean or per leaf segment)
            @param VPD [kPa]                vapour pressure deficit
            @param ci_init [mol mol-1]      initial estimation of the leaf internal [CO2] (mean or per leaf segment)
            @param pl_init [cm]             initial estimation of leaf matric poential (mean or per leaf segment)
            @param soil_k [day-1]           optionally, soil conductivities can be prescribed per segment, 
                                            conductivity at the root surface will be limited by the value, i.e. kr = min(kr_root, k_soil)
            @param log                      indicates if omputed values should be printed in a text file  
            @return [cm] root xylem pressure per root system node  
            
            ATT: assumes that leaf temparature is constant and not influenced by gs
         """
        if log:
            logfile = open('solve_leuning.txt', "w")
            logfile.write('input parameters: ')
            logfile.write(repr(self.Param))

        stop = False; self.Qlight = Qlight; self.log = log
        self.VPD = VPD; self.Tl = Tl; self.sim_time = sim_time
        diff = 0; loop = 0
        self.TairC = Tl - 273.15; self.ci = ci_init
        self.Param['cs'] = cs

        self.initVals(N, sim_time)

        x_old = np.full(len(self.rs.nodes), 0)
        An_old = gco2_old = ci_old = pg_old = np.full(len(self.seg_indxs), 0.)
        outputFlux_old = np.full((len(self.rs.nodes) - 1), 0.)
        self.x = p_linit

        while(stop != True):
            loop += 1
            if log:
                logfile = open('solve_leuning.txt', "a")
                logfile.write("*" * 120)
                logfile.write('\nLoop n°' + repr(loop))
                logfile.close()

            if (isinstance(self.x, float) or isinstance(self.x, int)):  # first loop
                self.rxi = self.rxj = self.x
            else:
                self.rxi = np.array([self.x[s[0]] for s in self.seg_leaves])
                self.rxj = np.array([self.x[s[1]] for s in self.seg_leaves])

            self.calcAn()
            self.calcGs()  # also computes E and Jw
            self.calcCi()
            self.calcPsig()  # water potential of guard cell
            self.linearSystem(sim_time, sxx, cells, soil_k)  # C++ (see XylemFlux.cpp)
            Q = sparse.coo_matrix((np.array(self.aV), (np.array(self.aI), np.array(self.aJ))))
            Q = sparse.csr_matrix(Q)
            self.x = LA.spsolve(Q, self.aB, use_umfpack = False)  # no difference when upfpack set to false

            outputFlux = self.radial_fluxes(2.0, self.x, sxx, [], True);

            diffX = max(abs(np.divide((self.x - x_old), self.x, out = np.zeros_like(self.x), where = self.x != 0)))
            diffAn = max(abs(np.divide((self.An - An_old), self.An, out = np.zeros_like(self.An), where = self.An != 0)))  # (self.An-An_old)/self.An))

            diffGs = max(abs(np.divide((self.gco2 - gco2_old), gco2_old, out = np.zeros_like(gco2_old), where = gco2_old != 0)))  # (self.gco2-gco2_old)/self.gco2))
            diffCi = max(abs(np.divide((self.ci - ci_old), self.ci, out = np.zeros_like(self.ci), where = self.ci != 0)))  # (self.ci-ci_old)/self.ci))
            diffPsig = max(abs(np.divide((self.getPg - pg_old), self.getPg, out = np.zeros_like(self.getPg), where = self.getPg != 0)))  # (self.getPg-pg_old)/self.getPg))
            diffFlux = max(abs(np.divide((outputFlux - outputFlux_old), outputFlux, out = np.zeros_like(outputFlux), where = outputFlux != 0)))  # (outputFlux-outputFlux_old)/outputFlux))
            diffRel = np.array([diffX, diffAn, diffGs, diffCi, diffPsig, diffFlux])

            # diffRel = np.array([max(abs(diffX/self.x)),max(abs(diffAn/self.An)), max(abs(diffGs/self.gco2)),
             #           max(abs(diffCi/self.ci)), max(abs(diffPsig/self.getPg))])

            x_old = self.x
            An_old = self.An
            gco2_old = self.gco2
            ci_old = self.ci
            pg_old = self.getPg
            outputFlux_old = outputFlux
            if(verbose):
                print("max error: ", diff, "max error (rel): ", diffRel, ", loop n°", loop)
            if((loop > 1000) or (max(abs(diffRel)) < 0.001)):
                stop = True
            if log:
                logfile = open('solve_leuning.txt', "a")
                logfile.write('\nAn (mol CO2 m-2 s-1) matrix: \n' + repr(self.An))
                logfile.write('\nRd (mol CO2 m-2 s-1) matrix: \n' + repr(self.Rd))
                logfile.write('\nVc (mol CO2 m-2 s-1) matrix: \n' + repr(self.Vc))
                logfile.write('\nVj (mol CO2 m-2 s-1) matrix: \n' + repr(self.Vj))
                logfile.write('\ngco2 (mol CO2 m-2 s-1) matrix: \n' + repr(self.gco2))
                logfile.write('\nci (mol mol-1) matrix: \n' + repr(self.ci))
                logfile.write('\nci/cs (-) matrix: \n' + repr([number / self.Param['cs'] for number in self.ci]))
                logfile.write('\np_l (cm) matrix: \n' + repr(self.x[self.leaf_nodes]))
                logfile.write('\nfw (-) matrix' + repr(self.fw) + '\n')
                logfile.write('\np_g (cm) matrix: \n' + repr(self.pg) + '\n')
                logfile.write('\nJw (cm3 cm-2 d-1) matrix: \n' + repr(self.Jw) + '\n')
                # logfile.write('\nE (cm3 d-1) matrix'+ repr(self.E)+'\n')
                logfile.close()
        res = max(abs(diffRel))
        if(verbose):
            print('leuning computation module stopped after {} trials.'.format(loop))
            print('Maximum absolute difference calculated at the last two trials for')
            print('psi:', diffX, 'An:', diffAn, 'gco2:', diffGs, 'ci:', diffCi, "guard cell's matric potential", diffPsig)
            print('max relative error in a cell:{}'.format(res))
        if log:
            logfile = open('solve_leuning.txt', "a")
            logfile.write('leuning computation module stopped after {} trials. Sum of absolute difference between leaf matric potential calculated at the last two trials: {}'.format(loop, diff))
            logfile.close()
            print('you can find the simulation log in the file solve_leuning.txt')
        return self.x

    def initVals(self, N, age):
        self.initStructandN(N, age)
        self.initVcVjRd()

    def initStructandN(self, N, age):
        self.leaf_nodes = self.get_nodes_index(4)
        self.segments = self.get_segments()
        self.seg_radii = np.array(self.rs.radii)
        self.organTypes = np.array(self.rs.organTypes)
        self.seg_length = np.array(self.rs.segLength())
        self.seg_indxs = np.array(self.get_segments_index(4))
        radii_leaf = self.seg_radii[self.seg_indxs]
        length_leaf = self.seg_length[self.seg_indxs]
        mass = (np.pi * pow(radii_leaf, 2) * length_leaf) * self.Param['rho_p']  # in g
        totArea = 0.0001 * (2 * np.pi * radii_leaf * length_leaf + 2 * np.pi * pow(radii_leaf, 2))  # m2

        if N == []:
            self.getNc()  # module to get N concentration in leaf. TODO: use N flow instead
        else:
            self.Nc = N  # in % or g N per g DM

        Np = ((self.Nc / 100) * mass / totArea) / 14  # in mol/m2
        self.Chl = self.Param['Chl/N'] * Np * 893.51 / 10  # chlorophyll a content in microg cm-2, evans1989
        self.seg_leaves = self.get_segments()[self.seg_indxs]
        self.sideArea = (2 * np.pi * radii_leaf * length_leaf)  # cm2
        st = 2
        kr = np.array([self.kr_f(age, st, 4, leafn) for leafn in range(len(self.seg_indxs))])
        kx = self.kx_f(age, st, 4)
        a = radii_leaf
        self.l = length_leaf
        self.f = -2 * a * np.pi * kr
        self.tau = np.sqrt(2 * a * np.pi * kr / kx)
        self.d = np.exp(-self.tau * self.l) - np.exp(self.tau * self.l)

        self.es = 0.61078 * np.exp(17.27 * self.TairC / (self.TairC + 237.3))  # 2.338205 kPa

    def initVcVjRd(self):
        # carboxylation rate
        # #Vc25max
        Vcrefmax = (self.Param['VcmaxrefChl1'] * self.Chl + self.Param['VcmaxrefChl2']) * 1e-6  # mol m-2 s-1
        # #Vcmax
        expo1 = np.exp(self.Param['Eav'] / (self.Param['R'] * self.Param['Tref']) * (1 - self.Param['Tref'] / self.Tl))
        expo2 = np.exp((self.Param['S'] * self.Tl - self.Param['Edv']) / (self.Param['R'] * self.Tl))
        self.Vcmax = Vcrefmax * expo1 / (1 + expo2)  # Eq 11
        # #Vc
        self.Ko = self.Param['Ko_ref'] * np.exp(self.Param['Eao'] / (self.Param['R'] * self.Param['Tref']) * (1 - self.Param['Tref'] / self.Tl))  # Eq 9
        self.Kc = self.Param['Kc_ref'] * np.exp(self.Param['Eac'] / (self.Param['R'] * self.Param['Tref']) * (1 - self.Param['Tref'] / self.Tl))  # Eq 9
        self.delta = self.Param['gamma0'] * (1 + self.Param['gamma1'] * (self.Tl - self.Param['Tref']) + self.Param['gamma2'] * pow((self.Tl - self.Param['Tref']), 2))  # Eq 10

        # electron transport rate
        # #Jrefmax
        Jrefmax = Vcrefmax * self.Param['a3']  # Eq 25
        # #Jmax
        expo1 = np.exp(self.Param['Eaj'] / (self.Param['R'] * self.Param['Tref']) * (1 - self.Param['Tref'] / self.Tl))
        expo2 = np.exp((self.Param['S'] * self.Tl - self.Param['Edj']) / (self.Param['R'] * self.Tl))
        Jmax = np.minimum(Jrefmax * expo1 / (1 + expo2), Jrefmax)  # Eq 24
        # #J
        coefa = self.Param['theta']
        coefb = -(self.Param['alpha'] * self.Qlight + Jmax)
        coefc = self.Param['alpha'] * self.Qlight * Jmax
        dis = pow(coefb, 2) - (4 * coefa * coefc)
        self.J = ((-coefb - np.sqrt(dis)) / (2 * coefa))

        self.Rd = self.Param['Rd_ref'] * np.exp(self.Param['Eard'] / (self.Param['R'] * self.Param['Tref']) * (1 - self.Param['Tref'] / self.Tl))
        if self.log:
            logfile = open('solve_leuning.txt', "a")
            logfile.write(' Chl ' + repr(self.Chl) + ' Vcrefmax ' + repr(Vcrefmax) + ', Vcmax ' + repr(self.Vcmax) +
            ', Ko ' + repr(self.Ko) + ', Kc ' + repr(self.Kc) + ', Γ* ' + repr(self.delta) +
            ', Jrefmax ' + repr(Jrefmax) + ', Jmax ' + repr(Jmax) + ', J ' + repr(self.J))
            logfile.close()

    def getNc(self):
        """ give N concentration in leaf. TODO: use N flow in plant instead
        """
        seg_aboveground = self.segments[np.where(self.organTypes > 2)  ]
        radii_aboveground = self.seg_radii[np.where(self.organTypes > 2)]
        length_aboveground = self.seg_length[np.where(self.organTypes > 2)]
        DM_aboveground = 0
        get_y_node = lambda vec: np.array(vec)[1]
        get_x_node = lambda vec: np.array(vec)[0]
        nodes_aboveground = np.concatenate((self.get_nodes_organ_type(3), self.get_nodes_organ_type(4)))
        nodesy = np.array([get_y_node(xi) for xi in self.segments])
        nodesx = np.array([get_x_node(xi) for xi in self.segments])
        surface = (max(nodesx) - min(nodesx)) * (max(nodesy) - min(nodesy)) * 1e-8  # cm2 -> ha
        DM_aboveground = sum((np.pi * np.power(radii_aboveground, 2) * length_aboveground) * self.Param['rho_p'] * 1e-6)  # g -> tonne,  #Eq 19
        DM_aboveground = DM_aboveground / surface  # Eq 19
        Nc = (DM_aboveground <= self.Param['DMcrit']) * 4.4 + (~(DM_aboveground <= self.Param['DMcrit'])) * (self.Param['anc'] * pow(DM_aboveground, -self.Param['bnc']))  # Eq 18
        self.Nc = Nc

    def calcAn(self):
        """
            fills the net assimilation and dark respiration vectors
            @param Q_input[mol photons m-2 s-1]     absorbed photon irradiance (mean or per leaf segment)
            @param Tl_inpu [K]                      leaf temperature (mean or per leaf segment)
            @param log                              indicates if omputed values should be printed in a text file  
            size = number of leaf segments
            vector is ordered like the leaf segments
        """
        # carboxylation rate
        self.Vc = np.minimum(np.maximum(self.Vcmax * (self.ci - self.delta) / (self.ci + self.Kc * (1 + self.Param['oi'] / self.Ko)), 0), self.Vcmax)  # Eq 8
        # electron transport rate
        # #Vj
        eps = 0
        if not isinstance(self.ci, float):
            eps = np.full(len(self.ci), 0.)
            eps[np.where([self.ci == 2 * self.delta])[0]] = 0.001 * self.delta  # to avoid a division by 0
        elif self.ci == 2 * self.delta:
            eps = 0.001 * self.delta
        self.Vj = np.maximum(self.J / 4 * (self.ci - self.delta) / (self.ci - 2 * self.delta + eps), 0)  # Eq 22
        # An and delta_gco2
        self.An = np.minimum(self.Vc, self.Vj) - self.Rd  # Eq 6
        self.deltagco2 = (self.delta + self.Kc * self.Rd * (1 + self.Param['oi'] / self.Ko) / self.Vcmax) / (1 - self.Rd / self.Vcmax)

    def calcGs(self):
        """ fills the stomatal conductance vector
            @param VPD_input [kPa]          vapour pressure deficit
            @param p_l_input [cm]           leaf matric poential (mean or per leaf segment)
            @param log                      indicates if omputed values should be printed in a text file  
            size = number of leaf segments
            vector is ordered like the leaf segments
        """
        p_lMPa = (self.rxi + self.rxj) * 0.5 * 0.0000978  # cm => MPa
        self.fw = self.Param['fwr'] + (1 - self.Param['fwr']) * np.exp(-np.exp(-self.Param['sh'] * (p_lMPa - self.Param['p_lcrit']) * 10228))  # Eq 5
        self.gco2 = self.Param['g0'] + self.fw * self.Param['a1'] * (self.An + self.Rd) / (self.ci - self.deltagco2)  # tuzet2003
        HRleaf = np.exp(self.Param['Mh2o'] * p_lMPa / (self.Param['rho_h2o'] * self.Param['R'] * self.Tl))  # fractional relative humidity in the intercellular spaces
        ea = self.es - self.VPD
        ea_leaf = self.es * HRleaf
        # VPDL = ea_leaf - ea
        self.Jw = (self.gco2 * 1.6) * (ea_leaf - ea) / self.Param['Patm'] * self.Param['Mh2o'] / self.Param['rho_h2o'] * 24 * 3600e-4  # in cm3 cm-2 d-1
        # self.E = self.Area * self.Jw

    def calcCi(self):
        """ give leaf segments internal [CO2] in mol mol-1 as matrix
        """
        self.ci = (self.Param['cs'] * self.Param['a1'] * self.fw + self.deltagco2) / (1 + self.Param['a1'] * self.fw)  # Eq 26

    def calcPsig(self):
        """ 
            fills the guard cell water potential vector
            @param Q_input[mol photons m-2 s-1]     absorbed photon irradiance (mean or per leaf segment) 
            size = number of leaf segments
            vector is ordered like the leaf segments
        """
        self.pg = -((self.sideArea * self.Jw) / (-self.f * (1. / (self.tau * self.d)) * (2. - np.exp(-self.tau * self.l) - np.exp(self.tau * self.l))) - (self.rxi + self.rxj)) / 2

    @property
    def getPg(self):
        return np.array(self.pg)
