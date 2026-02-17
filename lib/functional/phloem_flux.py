import timeit


import numpy as np
import json

import plantbox as pb
from plantbox import PhloemFlux
from functional.Photosynthesis import PhotosynthesisPython 



class PhloemFluxPython(PhloemFlux, PhotosynthesisPython):
    """  wrapper for photosynthesis
       
    """

    def __init__(self, plant_, params, psiXylInit, ciInit):
        """ @param mp is a pb.MappedPlant
            @param params are the hydraulic parameters
            @param psiXylInit [cm] is the initial guess of plant water potential [cm] for the fixed point iteration
            @param ciInit [mol mol-1] is the initial guess of leaf air CO2 partial pressure [-] for the fixed point iteration            
        """
        PhloemFlux.__init__( self,plant_,params, psiXylInit, ciInit)
        PhotosynthesisPython.__init__( self,plant_, params, psiXylInit, ciInit)
        self.reset()
        # self.update_outputs()
        
    def reset(self): # TODO: check
        self.Q_Rm      = np.array([])
        self.Q_Gr      = np.array([])
        self.Q_Exud    = np.array([])
        self.Q_ST      = np.array([])
        self.C_ST_np   = np.array([])    
        self.Q_meso    = np.array([])   
        self.C_meso    = np.array([])  
        self.Q_out = self.Q_out * 0
        self.Nt   = len(self.plant.nodes)
        self.Q_Rmbu      = np.array([])
        self.Q_Grbu      = np.array([])
        self.Q_Exudbu    = np.array([])
        self.Q_STbu    = np.array([])
        self.Q_mesobu    = np.array([])
        self.Ntbu = 0
        self.Q_ST_i        = np.array([])
        self.Q_Rm_i        = np.array([])
        self.Q_Gr_i        = np.array([])
        self.Q_Exud_i      = np.array([])
                
    def solve_phloem_flow(self, simDuration, dt, TairC, verbose = False, outputfile = "outputs.txt" ):
        self.startPM(simDuration, simDuration + dt, 1, ( TairC + 273.15) , verbose, outputfile )
        self.Nt = len(self.plant.nodes)        
        Q_out = np.array(self.Q_out) * 1e-3 * 12 # mmol Suc => mol C
        self.Q_ST    = np.array(Q_out[0:self.Nt])          #sieve tube sucrose content
        self.Q_meso  = np.array(Q_out[self.Nt:(self.Nt*2)])     #mesophyll sucrose content
        self.Q_Rm    = np.array(Q_out[(self.Nt*2):(self.Nt*3)]) #sucrose used for maintenance respiration
        self.Q_Exud  = np.array(Q_out[(self.Nt*3):(self.Nt*4)]) #sucrose used for exudation
        self.Q_Gr    = np.array(Q_out[(self.Nt*4):(self.Nt*5)]) #sucrose used for growth and growth respiration
        
        #self.Ntbu = len(self.Q_STbu)
        self.Q_STbu       =   np.concatenate((self.Q_STbu, np.full(self.Nt - self.Ntbu, 0.)))
        self.Q_Rmbu       =   np.concatenate((self.Q_Rmbu, np.full(self.Nt - self.Ntbu, 0.)))
        self.Q_Grbu       =   np.concatenate((self.Q_Grbu, np.full(self.Nt - self.Ntbu, 0.))) 
        self.Q_Exudbu     =   np.concatenate((self.Q_Exudbu, np.full(self.Nt - self.Ntbu, 0.))) 
            
        self.Q_ST_i        = self.Q_ST      - self.Q_STbu #in the sieve tubes
        self.Q_Rm_i        = self.Q_Rm      - self.Q_Rmbu #for maintenance
        self.Q_Gr_i        = self.Q_Gr      - self.Q_Grbu #for growth
        self.Q_Exud_i      = self.Q_Exud    - self.Q_Exudbu #for exudation
        #self.Q_out_i       = self.Q_Rm_i    + self.Q_Exud_i      + self.Q_Gr_i #total usage
                    
        volST   = np.array(self.vol_ST)         #sieve tube volume
        volMeso   = np.array(self.vol_Meso)      #mesophyll volume  
        self.C_ST_np    = np.array(self.C_ST)    
        self.C_meso  = self.Q_meso/volMeso  
        
        self.Ntbu = self.Nt
        self.Q_STbu       =   self.Q_ST.copy()
        self.Q_Rmbu       =   self.Q_Rm.copy()
        self.Q_Grbu       =   self.Q_Gr.copy() 
        self.Q_Exudbu     =   self.Q_Exud.copy()
    
    def update_outputs(self):    
        self.outputs_options = {
                "sieve tube concentration":self.C_ST_np, 
                "sieve tube content":self.Q_ST, 
                "mesophyll concentration":self.C_meso, 
                "mesophyll content":self.Q_meso, 
                "maintenance respiration":self.Q_Rm,
                "exudation":self.Q_Exud, 
                "growth":self.Q_Gr}
        self.outputs_options_last = {
                "sieve tube concentration":self.C_ST_np, 
                "sieve tube content":self.Q_ST, 
                "mesophyll concentration":self.C_meso, 
                "mesophyll content":self.Q_meso,
                "maintenance respiration":self.Q_Rm_i,
                "exudation":self.Q_Exud_i, 
                "growth":self.Q_Gr_i}
        
    def get_phloem_data_list(self): # TODO: complete
        return self.outputs_options.keys()
        
    def get_phloem_data(self, data, last = False, doSum = False):
        self.update_outputs()
        if last:
            outputs = self.outputs_options_last[data]
        else:
            outputs = self.outputs_options[data]
        if doSum:
            outputs = sum(outputs)
        return outputs
        
    
    def getPsiAir(self,RH, TairC):#constants are within photosynthesys.h
        return np.log(RH) * self.rho_h2o * self.R_ph * (TairC + 237.3)/self.Mh2o * (1/0.9806806)  ; #in cm
     
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

    def get_suf(self, sim_time):
        """ calculates the surface uptake fraction [1] of the root system at simulation time @param sim_time [day]
            (suf is constant for age independent conductivities)  """
        segs = self.rs.segments
        nodes = self.rs.nodes
        p_s = np.zeros((len(segs),))
        for i, s in enumerate(segs):
            p_s[i] = -500 - 0.5 * (nodes[s.x].z + nodes[s.y].z)  # constant total potential (hydraulic equilibrium)
        rx = self.solve_neumann(sim_time, -1.e5, p_s, False)  # False: matric potential not given per cell (but per segment), high number to recuce spurious fluxes
        print("rx ", np.min(rx), np.max(rx), np.mean(rx))
        fluxes = self.segFluxes(sim_time, rx, p_s, approx = False, cells = False)  # cm3/day, simTime,  rx,  sx,  approx, cells
        # print("fluxes ", np.min(fluxes) / -1.e5, np.max(fluxes) / -1.e5, np.mean(fluxes) / -1.e5)
        return np.array(fluxes) / -1.e5  # [1]

    def get_mean_suf_depth(self, sim_time):
        """  mean depth [cm] of water uptake based suf """
        suf = self.get_suf(sim_time)
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

    def get_krs(self, sim_time, seg_ind = [0]):
        """ calculatets root system conductivity [cm2/day] at simulation time @param sim_time [day] 
        if there is no single collar segment at index 0, pass indices using @param seg_ind, see find_base_segments        
        """
        segs = self.rs.segments
        nodes = self.rs.nodes
        p_s = np.zeros((len(segs),))
        for i, s in enumerate(segs):
            p_s[i] = -500 - 0.5 * (nodes[s.x].z + nodes[s.y].z)  # constant total potential (hydraulic equilibrium)
        rx = self.solve_dirichlet(sim_time, -15000, 0., p_s, cells = False)
        jc = 0
        for i in seg_ind:
            jc -= self.axial_flux(i, sim_time, rx, p_s, [], cells = False, ij = True)
        krs = jc / (-500 - 0.5 * (nodes[segs[0].x].z + nodes[segs[0].y].z) - rx[0])
        return krs , jc

    def get_eswp(self, sim_time, p_s):
        """ calculates the equivalent soil water potential [cm] at simulation time @param sim_time [day] for 
        the soil matric potential @param p_s [cm] given per cell """
        segs = self.rs.segments
        nodes = self.rs.nodes
        seg2cell = self.rs.seg2cell
        suf = self.get_suf(sim_time)
        eswp = 0.
        for i, s in enumerate(segs):
            eswp += suf[i] * (p_s[seg2cell[i]] + 0.5 * (nodes[s.x].z + nodes[s.y].z))  # matric potential to total potential
        return eswp

    def kr_f(self, age, st, ot = 2 , numleaf = 2, seg_ind = 0):
        """ root radial conductivity [1 day-1] for backwards compatibility """
        return self.kr_f_cpp(seg_ind, age, st, ot, numleaf)  # kr_f_cpp is XylemFlux::kr_f

    def kx_f(self, age, st, ot = 2, seg_ind = 0):
        """ root axial conductivity [cm3 day-1]  for backwards compatibility """
        return self.kx_f_cpp(seg_ind, age, st, ot)  # kx_f_cpp is XylemFlux::kx_f
            
    def write_phloem_parameters(self, filename="phloem_parameters"):
        """Write phloem flow module parameters to a JSON file."""
        parameters = {
            "InitialValues": {
                "initValST": {"value": self.initValST, "description": "Initial concentration in sieve tube"},
                "initValMeso": {"value": self.initValMeso, "description": "Initial concentration in mesophyll"},
                "withInitVal": {"value": self.withInitVal, "description": "Use initial values"}
            },
            "Growth": {
                "psi_osmo_proto": {"value": self.psi_osmo_proto, "unit": "cm", "description": "Osmotic potential in protophloem"},
                "psiMin": {"value": self.psiMin, "unit": "cm", "description": "Minimum water potential for growth"},
                "leafGrowthZone": {"value": self.leafGrowthZone, "unit": "cm", "description": "Leaf growth zone length"},
                "Gr_Y": {"value": self.Gr_Y, "description": "Growth efficiency"},
                "StemGrowthPerPhytomer": {"value": self.StemGrowthPerPhytomer, "description": "Growth per phytomer"},
                "useCWGr": {"value": self.useCWGr, "description": "Use C and W limited growth"}
            },
            "SieveTube": {
                "Vmaxloading": {"value": self.Vmaxloading, "unit": "mmol cm-1 d-1", "description": "Max sucrose loading"},
                "CSTimin": {"value": self.CSTimin, "description": "Minimum sucrose threshold"},
                "beta_loading": {"value": self.beta_loading, "description": "Feedback effect of C_ST"},
                "Mloading": {"value": self.Mloading, "description": "Michaelis-Menten coefficient for loading"},
                "C_targ": {"value": self.C_targ, "unit": "mmol Suc cm-3", "description": "Sucrose target concentration"},
                "Q10": {"value": self.Q10, "description": "Q10 value for respiration"},
                "TrefQ10": {"value": self.TrefQ10, "unit": "°C", "description": "Reference temperature for Q10"},
                "KMfu": {"value": self.KMfu, "description": "Michaelis-Menten coefficient for sucrose usage"},
                "k_mucil": {"value": self.k_mucil, "unit": "d-1", "description": "Decay rate of mucilage"},
                "k_mucil_": {"value": self.k_mucil_, "description": "Vector of mucilage decay rates"},
                "Vmax_S_ST": {"value": self.Vmax_S_ST, "unit": "mmol Suc d-1 cm-3", "description": "Max sucrose usage"},
                "kM_S_ST": {"value": self.kM_S_ST, "unit": "mmol Suc cm-3", "description": "Michaelis-Menten constant for ST"},
                "kHyd_S_ST": {"value": self.kHyd_S_ST, "unit": "d-1", "description": "Sucrose hydrolysis rate"},
                "k_S_ST": {"value": self.k_S_ST, "unit": "d-1", "description": "Sucrose loss rate"},
                "update_viscosity_": {"value": self.update_viscosity, "description": "Update viscosity"},
                "usePsiXyl": {"value": self.usePsiXyl, "description": "Use xylem water potential to get total phloem potential"}
            },
            "Mesophyll": {
                "C_targMesophyll": {"value": self.C_targMesophyll, "unit": "mmol Suc cm-3", "description": "Target sucrose concentration"},
                "Vmax_S_Mesophyll": {"value": self.Vmax_S_Mesophyll, "unit": "mmol Suc d-1 cm-3", "description": "Max sucrose usage in mesophyll"},
                "kM_S_Mesophyll": {"value": self.kM_S_Mesophyll, "unit": "mmol Suc cm-3", "description": "Michaelis-Menten constant"},
                "kHyd_S_Mesophyll": {"value": self.kHyd_S_Mesophyll, "unit": "d-1", "description": "Hydrolysis rate"},
                "k_S_Mesophyll": {"value": self.k_S_Mesophyll, "unit": "d-1", "description": "Loss rate"},
                "surfMeso": {"value": self.surfMeso, "unit": "cm2", "description": "Cross-sectional area of mesophyll"},
                "sameVolume_meso_seg": {"value": self.sameVolume_meso_seg, "description": "Same volume for mesophyll and segment"},
                "sameVolume_meso_st": {"value": self.sameVolume_meso_st, "description": "Same volume for sieve tube and mesophyll"},
            },
            "PerType": { # TODO: move unit change to wrapper and use mal and cm in python upper layer
                "Across_st": {"value": self.Across_st, "unit": "cm2", "description": "effect of the sucrose content on maintenance respiration"},
                "kr_st": {"value": self.kr_st, "unit": "mmol hPa-1 day-1", "description": "effect of the sucrose content on maintenance respiration"},
                "kx_st": {"value": self.kx_st, "unit": "cm3 hPa-1 day-1", "description": "effect of the sucrose content on maintenance respiration"},
                "Krm2": {"value": self.krm2v, "unit": "-", "description": "effect of the sucrose content on maintenance respiration"},
                "Krm1": {"value": self.krm1v, "unit": "-", "description": "effect of the sucrose content on maintenance respiration"},
                "Rho_s": {"value": self.rhoSucrose, "unit": "mmol Suc cm-3", "description": "sucrose density per organ type"},
                "Rmax_st": {"value": self.Rmax_st, "unit": "cm d-1", "description": "maximum growth rate when water and carbon limitation is activated"}
            },
            "Soil": {
                "DefaultC": {"value": self.CsoilDefault, "unit": "mmol Suc cm-3", "description": "dummy value for soil concentration"}
            },
            "Solver": {
                "atol": {"value": self.atol, "description": "Absolute tolerance"},
                "rtol": {"value": self.rtol, "description": "Relative tolerance"},
                # "solver": {"value": self.solver, "description": "Solver type"},
                "doTroubleshooting": {"value": self.doTroubleshooting, "description": "Enable troubleshooting"}
            }
        }

        with open(filename + ".json", "w+") as f:
            json.dump(parameters, f)

    def read_phloem_parameters(self, filename):
        """Read phloem flow module parameters from a JSON file."""
        with open(filename + ".json", "r") as f:
            parameters = json.load(f)

        self.initValST = parameters["InitialValues"]["initValST"]["value"]
        self.initValMeso = parameters["InitialValues"]["initValMeso"]["value"]
        self.withInitVal = parameters["InitialValues"]["withInitVal"]["value"]

        self.psi_osmo_proto = parameters["Growth"]["psi_osmo_proto"]["value"]
        self.psiMin = parameters["Growth"]["psiMin"]["value"]
        self.leafGrowthZone = parameters["Growth"]["leafGrowthZone"]["value"]
        self.Gr_Y = parameters["Growth"]["Gr_Y"]["value"]
        self.StemGrowthPerPhytomer = parameters["Growth"]["StemGrowthPerPhytomer"]["value"]
        self.useCWGr = parameters["Growth"]["useCWGr"]["value"]

        self.Vmaxloading = parameters["SieveTube"]["Vmaxloading"]["value"]
        self.CSTimin = parameters["SieveTube"]["CSTimin"]["value"]
        self.beta_loading = parameters["SieveTube"]["beta_loading"]["value"]
        self.Mloading = parameters["SieveTube"]["Mloading"]["value"]
        self.C_targ = parameters["SieveTube"]["C_targ"]["value"]
        self.Q10 = parameters["SieveTube"]["Q10"]["value"]
        self.TrefQ10 = parameters["SieveTube"]["TrefQ10"]["value"]
        self.KMfu = parameters["SieveTube"]["KMfu"]["value"]
        self.k_mucil = parameters["SieveTube"]["k_mucil"]["value"]
        self.k_mucil_ = parameters["SieveTube"]["k_mucil_"]["value"]
        self.Vmax_S_ST = parameters["SieveTube"]["Vmax_S_ST"]["value"]
        self.kM_S_ST = parameters["SieveTube"]["kM_S_ST"]["value"]
        self.kHyd_S_ST = parameters["SieveTube"]["kHyd_S_ST"]["value"]
        self.k_S_ST = parameters["SieveTube"]["k_S_ST"]["value"]
        self.update_viscosity_ = parameters["SieveTube"]["update_viscosity_"]["value"]
        self.usePsiXyl = parameters["SieveTube"]["usePsiXyl"]["value"]

        self.C_targMesophyll = parameters["Mesophyll"]["C_targMesophyll"]["value"]
        self.Vmax_S_Mesophyll = parameters["Mesophyll"]["Vmax_S_Mesophyll"]["value"]
        self.kM_S_Mesophyll = parameters["Mesophyll"]["kM_S_Mesophyll"]["value"]
        self.kHyd_S_Mesophyll = parameters["Mesophyll"]["kHyd_S_Mesophyll"]["value"]
        self.k_S_Mesophyll = parameters["Mesophyll"]["k_S_Mesophyll"]["value"]
        self.surfMeso = parameters["Mesophyll"]["surfMeso"]["value"]
        self.sameVolume_meso_seg = parameters["Mesophyll"]["sameVolume_meso_seg"]["value"]
        self.sameVolume_meso_st = parameters["Mesophyll"]["sameVolume_meso_st"]["value"]
        
        self.setKrm2(parameters["PerType"]["Krm2"]["value"]) 
        self.setKrm1(parameters["PerType"]["Krm1"]["value"]) 
        self.setRhoSucrose(parameters["PerType"]["Rho_s"]["value"]) 
        self.setRmax_st(parameters["PerType"]["Rmax_st"]["value"]) 
        self.setKr_st(parameters["PerType"]["kr_st"]["value"]) 
        self.setKx_st(parameters["PerType"]["kx_st"]["value"]) 
        self.setAcross_st(parameters["PerType"]["Across_st"]["value"]) 
                    
        self.CsoilDefault = parameters["Soil"]["DefaultC"]["value"]

        self.atol = parameters["Solver"]["atol"]["value"]
        self.rtol = parameters["Solver"]["rtol"]["value"]
        # self.solver = parameters["Solver"]["solver"]["value"]
        self.doTroubleshooting = parameters["Solver"]["doTroubleshooting"]["value"]

    @staticmethod
    def convert_(x, dtype = np.float64):
        """ not used anymore (?) """
        return np.array(list(map(lambda x: np.array(x, dtype), x)), dtype)  # is there a better way?
