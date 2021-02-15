import timeit
import math

import numpy as np
from scipy import sparse
import scipy.sparse.linalg as LA
import collections

import plantbox as pb
from plantbox import XylemFlux
import rsml_reader as rsml  # todo
from xylem_flux import XylemFluxPython  # Python hybrid solver

class Leuning(XylemFluxPython):
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
        self.An = [] #
        self.Rd = []#
        self.ci = []#
        self.gco2 =  []# stomatal aperture for CO2
        self.Nc = 0 # in % or g N per g DM
        self.Param = { #list of parameters in alphabetical order
        'a1' :6, #(-) fitting parameter for gco2
        'a2': 1.6,#(-)
        'a3':2,
        'alpha': 0.2,
        'anc': 5.35,
        'arn': 1.8351,
        'bnc': 0.44,
        'brn' : -1.5075,
        'c': 7.72*5, #added *5 with regards to Lu2020
        'cs':350e-6, #example from Dewar2002
        'd': 0.04,
        'D0': 1.67, #kPa
        'DMcrit': 1.55,
        'Eac': 59430,
        'Eaj': 37000,
        'Eao': 36000,
        'Eav': 58520,
        'Edj':  220000,
        'Edv': 220000,
        'fwr': 0,
        'gamma0': 28e-6,
        'gamma1': 0.0509,
        'gamma2': 0.001,
        'Ko_ref': 256e-3,
        'Kc_ref': 302e-6,
        'Mh2o': 18,
        'Ncmin': 4.4,
        'oi': 210e-3,
        'p_lcrit': -5500,
        'R': 8.31,
        'rho_h2o': 1,
        'rho_p': 0.35,
        'S': 700,
        'sh': 1e100,
        'Tref': 293.2,
        'theta': 0.9,
        }
    
    def solve_leuning(self, sim_time :float,sxx, cells :bool, Qlight ,VPD: float,Tl, p_linit,ci_init,  soil_k = [], log = True) :
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
         """
        if log:
            logfile = open('../../../tutorial/examples/python/results/solve_leuning.txt', "w")
            logfile.write('input parameters: ')
            logfile.write(repr(self.Param))
        
        organTypes = self.get_organ_types()
        leaf_nodes = self.get_nodes_index(4)
        stop = False
        diff = 0
        loop = 0
        p_l = p_linit 
        self.gs = []
        self.gco2 = []
        self.An = []
        self.Rd = [] # reset the arrays
        self.ci = ci_init
        Indx = np.concatenate((self.get_organ_nodes_tips()))
        value = np.full(len(Indx), 0) #set 0 axial fux at tip of stems, roots and leaves
        x_old = np.full(len(self.rs.nodes), 0)
        while(stop != True):   
            loop +=1 
            if log:
                logfile.write("*" * 120)
                logfile.write('\nLoop n°'+repr(loop))
                logfile.close()
            self.getNc() #module to get N concentration in leaf. TODO: use N flow instead
            self.calcAn(Qlight, Tl, log) 
            self.calcGs(VPD,  p_l, log) 
            self.calcCi() 
            self.linearSystem(sim_time, sxx, cells, soil_k) # C++ (see XylemFlux.cpp)
            Q = sparse.coo_matrix((np.array(self.aV), (np.array(self.aI), np.array(self.aJ))))
            Q = sparse.csr_matrix(Q)
            Q, b = self.bc_neumann(Q, self.aB, Indx, value)  # cm3 day-1
            x = LA.spsolve(Q, b, use_umfpack = True)  # direct
            p_l = x[leaf_nodes]
            diff = sum(np.sqrt((x-x_old)**2))
            x_old = x
            if(loop > 1000 or diff < 1.e-5):
                stop = True
            if log:
                logfile = open('../../../tutorial/examples/python/results/solve_leuning.txt', "a")
                logfile.write('\nAn (mol CO2 m-2 s-1) matrix: \n' +repr(self.An))
                logfile.write('\nRd (mol CO2 m-2 s-1) matrix: \n' +repr(self.Rd))
                logfile.write('\ngco2 (mol CO2 m-2 s-1) matrix: \n' +repr(self.gco2))
                logfile.write('\nci (mol mol-1) matrix: \n' +repr(self.ci))
                logfile.write('\nci/cs (-) matrix: \n' +repr([number / self.Param['cs'] for number in self.ci]))
                logfile.write('\np_l (cm) matrix: \n' +repr(p_l))
                logfile.write('\nkr (d-1) matrix: \n' +repr(self.gs)+'\n')
        print('leuning computation module stopped after {} trials. Sum of absolute difference between leaf matric potential calculated at the last two trials: {}'.format(loop, diff))
        if log:
            logfile.write('leuning computation module stopped after {} trials. Sum of absolute difference between leaf matric potential calculated at the last two trials: {}'.format(loop, diff))
            print('you can find the simulation log in the file ../../../tutorial/examples/python/results/solve_leuning.txt')
        return x
        
    def calcAn( self,Q_input, Tl_input, log ):
        """
            fills the net assimilation and dark respiration vectors
            @param Q_input[mol photons m-2 s-1]     absorbed photon irradiance (mean or per leaf segment)
            @param Tl_inpu [K]                      leaf temperature (mean or per leaf segment)
            @param log                              indicates if omputed values should be printed in a text file  
            size = number of leaf segments
            vector is ordered like the leaf segments
        """
        self.An = [] # reset the array
        self.Rd = []
        seg_radii = np.array(self.rs.radii)
        organTypes = np.array(self.rs.organTypes)
        seg_length = np.array(self.segLength())
        radii_leaf = seg_radii[np.where(organTypes == 4)]  
        length_leaf = seg_length[np.where(organTypes == 4)]
        for si in range(len(length_leaf)) :
            if(isinstance(Q_input,(float,int))):
                Q = Q_input
            else:
                Q = Q_input[si]
            if(isinstance(Tl_input,(float,int))):
                Tl = Tl_input
            else:
                Tl = Tl_input[si]
            if(isinstance(self.ci,(float,int))):
                ci = self.ci
            else:
                ci = self.ci[si]
            #carboxylation rate
            ##Vc25max
            mass = (math.pi * pow(radii_leaf[si],2) * length_leaf[si])* self.Param['rho_p'] #in g
            Area = 0.0001 * (2 * math.pi * radii_leaf[si]* length_leaf[si] + 2 * math.pi * pow(radii_leaf[si] , 2)) 
            Np = (self.Nc/100) * mass / Area #in g/m2, Eq 21
            F_LNR = (self.Param['arn'] + self.Param['brn']/Np)*0.16 #Eq 16
            Rub = F_LNR*Np*6.25 #g/m² #Eq 15
            Rub_mol = Rub * (1e6/55000) # in micromol/m² #Eq 14
            K25cat = self.Param['c']/(1 + self.Param['d']*Rub_mol*8) #Eq 13
            Vcrefmax = K25cat * (8/550) * F_LNR * Np * 6.25 * 0.001  #Eq 12
            ##Vcmax
            expo1 = math.exp(self.Param['Eav'] /(self.Param['R']*self.Param['Tref'])*(1 - self.Param['Tref']/Tl))
            expo2 = math.exp((self.Param['S'] * Tl - self.Param['Edv'])/(self.Param['R'] * Tl))
            Vcmax = Vcrefmax * expo1 / (1 + expo2) #Eq 11
            ##Vc
            Ko = self.Param['Ko_ref'] * math.exp(self.Param['Eao']/(self.Param['R']*self.Param['Tref'])*(1-self.Param['Tref']/Tl)) #Eq 9
            Kc = self.Param['Kc_ref'] * math.exp(self.Param['Eac']/(self.Param['R']*self.Param['Tref'])*(1-self.Param['Tref']/Tl)) #Eq 9
            delta = self.Param['gamma0'] * (1+ self.Param['gamma1']*(Tl - self.Param['Tref']) + self.Param['gamma2']*pow((Tl - self.Param['Tref']),2) ) #Eq 10
            Vc = Vcmax * (ci - delta) / (ci + Kc*(1 + self.Param['oi']/Ko)) #Eq 8
            #electron transport rate
            ##Jrefmax
            Jrefmax = Vcrefmax * self.Param['a3'] #Eq 25
            ##Jmax
            expo1 = math.exp(self.Param['Eaj'] /(self.Param['R']*self.Param['Tref'])*(1 - self.Param['Tref']/Tl))
            expo2 = math.exp((self.Param['S'] * Tl - self.Param['Edj'])/(self.Param['R'] * Tl))
            Jmax = Jrefmax * expo1 / (1 + expo2) #Eq 24
            ##J
            coefa = self.Param['theta']
            coefb = -(self.Param['alpha'] * Q + Jmax)
            coefc = self.Param['alpha'] * Q * Jmax
            dis = pow(coefb,2) - (4*coefa*coefc)
            if (dis < 0):
                print('leuning::calcAn() : θJ^2-(αQ+J_max )J+αQJ_max=0, dis < 0, cannot solve the quadratic equation')
                raise
            if (coefa == 0):
                print('leuning::calcAn() : θJ^2-(αQ+J_max )J+αQJ_max=0, θ = 0, cannot solve the quadratic equation')
                raise
            J = ((-coefb+ math.sqrt(dis))/(2*coefa)) #Eq 23
            ##Vj
            Vj = J/4 * (ci - delta)/ (ci - 2 * delta) #Eq 22
            
            #An and Rd
            Rd = 0.015 * Vc #Eq 7
            An = max(min(Vc, Vj) - Rd, 0) #Eq 6
            self.An.append(An) #cm3/day
            self.Rd.append(Rd) #cm3/day
            if log:
                logfile = open('../../../tutorial/examples/python/results/solve_leuning.txt', "a")
                logfile.write('\nleaf seg n°' +repr(si))
                logfile.write(' Np '+ repr(Np)+ ', Vcrefmax '+ repr(Vcrefmax)+', Vcmax '+ repr(Vcmax)+
                ', Vc '+ repr(Vc)+', Ko '+ repr(Ko)+', Kc '+ repr(Kc)+ ', Γ* '+ repr(delta)+ 
                ', Jrefmax '+ repr(Jrefmax)+', Jmax '+ repr(Jmax)+', J '+ repr(J)+ ', Vj '+ repr(Vj))
                logfile.close()
            
    def calcGs( self,VPD_input, p_l_input, log):
        """ fills the stomatal conductance vector
            @param VPD_input [kPa]          vapour pressure deficit
            @param p_l_input [cm]           leaf matric poential (mean or per leaf segment)
            @param log                      indicates if omputed values should be printed in a text file  
            size = number of leaf segments
            vector is ordered like the leaf segments
        """
        gsleaf = []
        self.gs = []
        self.gco2 = []
        for si in range(len(self.get_segments_index(4))) :
            if(isinstance(VPD_input,(float,int))):
                VPD = VPD_input
            else:
                VPD = VPD_input[si]
            if(isinstance(p_l_input,(float,int))):
                p_l = p_l_input
            else:
                p_l = p_l_input[si]
            if(isinstance(self.ci,(float,int))):
                ci = self.ci
            else:
                ci = self.ci[si]
            fw = self.Param['fwr'] + (1- self.Param['fwr'])*math.exp(-math.exp(-self.Param['sh']*(p_l - self.Param['p_lcrit']))) #Eq 5
            gs =  fw * self.Param['a1'] * (self.An[si] + self.Rd[si]) / (ci*(1 + VPD/self.Param['D0'])) #Eq 1
            self.gco2.append(gs)
            Jwr = self.Param['a2'] *gs * self.Param['Mh2o'] / self.Param['rho_h2o'] * 3600e-4  #Eq 2 and 3
            kr = Jwr / abs(self.airPressure - p_l) #Eq 4
            gsleaf.append(kr)
            if log:
                logfile = open('../../../tutorial/examples/python/results/solve_leuning.txt', "a")
                logfile.write('\nleaf seg n°' +repr(si))
                logfile.write(' Jw (cm3 cm-2 d-1) '+ repr(Jwr))
                logfile.close()
        self.gs = gsleaf 
            
    def getNc(self):
        """ give N concentration in leaf. TODO: use N flow in plant instead
        """
        seg_radii = np.array(self.rs.radii)
        organTypes = np.array(self.rs.organTypes)
        seg_length = np.array(self.segLength())
        seg = np.array(self.get_segments())
        DM_aboveground = 0
        seg_aboveground = seg[np.where(organTypes > 2)  ]
        radii_aboveground = seg_radii[np.where(organTypes > 2)]  
        length_aboveground = seg_length[np.where(organTypes > 2)]  
        get_y_node = lambda vec : np.array(vec)[1]
        get_x_node = lambda vec : np.array(vec)[0]
        nodes_aboveground = np.concatenate((self.get_nodes_organ_type(3),self.get_nodes_organ_type(4)))
        nodesy = np.array([get_y_node(xi) for xi in seg])
        nodesx = np.array([get_x_node(xi) for xi in seg])
        surface = (max(nodesx) - min(nodesx)) * (max(nodesy) - min(nodesy)) * 1e-8 #cm2 -> ha
        for si in range(len(seg_aboveground)):
            DM_aboveground += (math.pi * pow(radii_aboveground[si],2) * length_aboveground[si])* self.Param['rho_p']*1e-6 # g -> tonne,  #Eq 19
        DM_aboveground = DM_aboveground/surface  #Eq 19
        Nc = (DM_aboveground <= self.Param['DMcrit'])*4.4 + (not(DM_aboveground <= self.Param['DMcrit'])) * (self.Param['anc']*pow(DM_aboveground,-self.Param['bnc']))  #Eq 18
        self.Nc = Nc
        
    def calcCi(self):
        """ give leaf segments internal [CO2] in mol mol-1 as matrix
        """
        self.ci = []
        for si in range(len(self.get_segments_index(4))) :
            ci = self.Param['cs'] - self.An[si]/self.gco2[si]  #Eq 26
            self.ci.append(ci)