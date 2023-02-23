import timeit
import math

import numpy as np
from scipy import sparse
import scipy.sparse.linalg as LA

import plantbox as pb
from plantbox import XylemFlux
import rsml_reader as rsml  # todo
from xylem_flux import XylemFluxPython  # Python hybrid solver

class StomataModel(XylemFluxPython):
    """  Hybrid flux solver (following Meunier et al.)
    
        C++ part is defined in CPlantBox XylemFlux.hh and XylemFlux.cpp
        
        Calculates water movement within the xylems, assuming a constant matric potential around each xylem segment,
        given for each segment, or for each soil cell. 
        
        The root surface flux is calculated exactly as in Meunier et al.        
    """

        
    def __init__(self, rs, PAR: float, VPD:float, TH:float, TL: float, Topt:float, psi1: float ,psi2: float, gmax:float):
        """ @param rs is either a pb.MappedPlant or a string containing a rsml filename
            sets the parameters for the calculation of gs and of the transpiration rate
            see lobet et al 2014 (planetMaiz )                          
        """
        super().__init__(rs)
        self.PAR = PAR; #[µmol/m2/sec]
        self.VPD = VPD;#MPa
        self.TH = TH;#*C
        self.TL = TL;#*C
        self.Topt = Topt;#*C
        self.psi1 = psi1;#MPa
        self.psi2 = psi2;#MPa
        self.gmax = gmax;#cm3/day
    
    def solve_neumann_gs(self, sim_time :float,sxx, cells :bool, PAR: float,VPD: float,Tair: float, p_linit,  soil_k = []) :
        """ solves the flux equations, with a neumann boundary condtion, see solve()
            @param sim_time [day]       needed for age dependent conductivities (age = sim_time - segment creation time)
            @param sxx [cm]             soil matric potentials given per segment or per soil cell
            @param cells                indicates if the matric potentials are given per cell (True) or by segments (False)
            @param soil_k [day-1]       optionally, soil conductivities can be prescribed per segment, 
                                        conductivity at the root surface will be limited by the value, i.e. kr = min(kr_root, k_soil)  
            @return [cm] root xylem pressure per root system node         
         """
        organTypes = self.get_organ_types()
        leaf_nodes = self.get_nodes_index(4)
        stop = False
        diff = 0
        loop = 0
        p_l = p_linit 
        Indx = np.concatenate((self.get_organ_nodes_tips()))
        value = np.full(len(Indx), 0) #set 0 axial fux at tip of stems, roots and leaves
        x_old = np.full(len(self.rs.nodes), 0)
        while(stop != True):    
            self.calcGs(PAR, VPD, Tair,  p_l) 
            print("gs", self.gs)
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
            loop +=1
        print('gs and rx computation module stopped after {} trials. Sum of absolute diference between rx calculated at the last two trials: {}'.format(loop, diff))
        return x
    
    def calcGs( self,PAR: float, VPD:float, Tair:float, p_leaf):
        """
            fills the stomatal conductance vector
            size = number of leaf segments
            vector is ordered like the leaf segments
        """
        gs = []
        for si in range(len(self.get_segments_index(4))) :
            if(isinstance(p_leaf,(float,int))):
                p_l = p_leaf*0.0000978
            else:
                p_l = p_leaf[si]*0.0000978 #1 cm of water = 0.0000978 MPa
            F_PAR = min(PAR/self.PAR,1)
            F_VPD = min(np.exp(-self.VPD*VPD),1)
            powTair = (self.TH - self.Topt)/(self.Topt - self.TL)
            divisorTair =((self.Topt - self.TL)*(self.TH - self.Topt))
            dividendTair =(Tair-self.TL)*(self.TH - Tair)
            F_Tair = (dividendTair / divisorTair) ** powTair
            F_pleaf = min(np.exp(-((-p_l/self.psi1)** self.psi2)),1)
            F_all = F_pleaf * F_Tair *F_VPD * F_PAR
            gs.append(self.gmax * F_all) #cm3/day
        self.gs = gs
            
        