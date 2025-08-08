import timeit


import numpy as np
import json
#import matplotlib.pyplot as plt

import plantbox as pb
from plantbox import Photosynthesis
from functional.PlantHydraulicModel import HydraulicModel_Meunier  

#import rsml_reader as rsml



class PhotosynthesisPython(Photosynthesis, HydraulicModel_Meunier):
    """  wrapper for photosynthesis
       
    """

    def __init__(self, mp, params, psiXylInit = -100, ciInit = 2e-4):
        """ @param mp is a pb.MappedPlant
            @param psiXylInit [cm] is the initial guess of plant water potential [cm] for the fixed point iteration
            @param ciInit [mol mol-1] is the initial guess of leaf air CO2 partial pressure [-] for the fixed point iteration            
        """
        Photosynthesis.__init__(self, mp, params, psiXylInit, ciInit)
        HydraulicModel_Meunier.__init__(self,mp, params)
        self.pCO2 = 350e-6
        
    def solve(self,sim_time:float, rsx:list,
                ea:float, es:float, PAR, TairC, cells = True, soil_k = [], 
                verbose = False, doLog = False, outputDirectory = "./results/"):
        """ Solves the stomatal regulation/photosynthesis + hydraulic model using the Meunier method 
            @param sim_time [day]   needed for age dependent conductivities (age = sim_time - segment creation time)
            @param rsx [cm]         soil matric potentials given per segment or per soil cell
            @param cells            indicates if the matric potentials are given per cell (True) or by segments (False)
            @param soil_k [day-1]   optionally, soil conductivities can be prescribed per segment, 
                                    conductivity at the root surface will be limited by the value, i.e. kr = min(kr_root, k_soil)  
            @param ea [kPa]         actual vapour pressure
            @param es [kPa]         reference vapour pressure
            @param PAR [mol photons cm-2 d-1]         absorbed photon irradiance per leaf segment
            @param TairC [°C]      leaf temperature (mean)
            @param verbose          print at runtime nothing (0), sparsly (1), many outputs (2)
            @param doLog            indicates if computed values should be printed in a text file (True) or not (False) 
            @param outputDirectory  where to save the log files
        """
        
        # send python data to cpp
        if isinstance(PAR, (int,float)):
            self.Qlight = [PAR / (24 * 3600) * 1e4] # to [mol photons m-2 s-1]
        elif isinstance(PAR, (type(np.array([])),type([]))):
            if isinstance(PAR, type([])):
                PAR = np.array(PAR)
            self.Qlight = PAR / (24 * 3600) * 1e4
        else:
            raise Exception(f'unexpected object type for PAR ({type(PAR)})')
            
        if isinstance(self.pCO2, (int,float)):
            self.cs = [self.pCO2] 
        elif isinstance(self.pCO2, (type(np.array([])),type([]))):
            self.cs = self.pCO2
        else:
            raise Exception(f'unexpected object type for pCO2 ({type(self.pCO2)})')
    
        if isinstance(TairC, (int,float)):
            TairK = [TairC + 273.15 ] 
        elif isinstance(TairC, (type(np.array([])),type([]))):
            if isinstance(TairC, type([])):
                TairC = np.array(TairC)
            TairK = TairC + 273.15
        else:
            raise Exception(f'unexpected object type for TairC ({type(TairC)})')
            
        self.solve_photosynthesis(sim_time = sim_time, sxx = rsx, cells = cells, 
                                        ea =  ea, es = es, TleafK =  TairK , soil_k = soil_k,  
                                verbose = verbose, doLog = doLog, outputDir = outputDirectory)
                                
    
    def get_leafBlade_Idx(self, ot = -1):
        """ return the global indexes of segments were the leaf blade area > 0
            useful, e.g., when looking at mean(An) 
            @param ot [-]   organ type of segments for which the elaf blade area is returned. 
                            does not  change the value of the array but its size
        """
        leafBlade = np.array(self.plant.leafBladeSurface)
        if ot != -1:
            organTypes = np.array(self.plant.organTypes)
            leafBlade = leafBlade[organTypes == ot]
        return np.where(leafBlade > 0)[0]
        
    def get_leafBlade_area(self, ot = -1):
        """ return the indexes of segments were the leaf blade area > 0 
            @param ot [-]   organ type of segments for which the elaf blade area is returned. 
                            does not  change the value of the array but its size
        """
        leafBlade = np.array(self.plant.leafBladeSurface)
        if ot == -1:
            return leafBlade
        else :
            organTypes = np.array(self.plant.organTypes)
            return np.array(leafBlade)[organTypes == ot]
        
    def get_transpiration(self):
        """ actual transpiration [cm3 day-1], calculated as the sum of all leaf radial fluxes"""
        return self.Ev

    def get_net_assimilation(self):
        """ net actual assimilation rate assimilation [mol CO2 d-1] """
        leafBlade = self.get_leafBlade_area(pb.leaf) * 2.
        An = self.get_net_assimilation_perleafBladeArea()
        return An * leafBlade
    
    def get_Vc(self):
        """ gross carboxilation-limited assimilation rate [mol CO2 d-1]"""
        leafBlade = self.get_leafBlade_area(pb.leaf) * 2.
        Vc = self.get_Vc_perleafBladeArea()
        return Vc * leafBlade
        
    def get_Vj(self):
        """ gross electron transport-limited assimilation rate [mol CO2 d-1] """
        leafBlade = self.get_leafBlade_area(pb.leaf) * 2.
        Vj = self.get_Vj_perleafBladeArea()
        return Vj * leafBlade
        
    def get_net_assimilation_perleafBladeArea(self):
        """ net actual assimilation rate per unit of surface [mol CO2 cm^-2 d-1] """
        return np.array(self.An) * 3600 * 24 / 1e4
        
    def get_Vc_perleafBladeArea(self):
        """ gross carboxilation-limited assimilation rate per unit of surface [mol CO2 cm^-2 d-1] """
        return np.array(self.Vc) * 3600 * 24 / 1e4
        
    def get_Vj_perleafBladeArea(self):
        """ gross electron transport-limited assimilation rate per unit of surface [mol CO2 cm^-2 d-1] """
        return np.array(self.Vj) * 3600 * 24 / 1e4
        
    def get_water_potential(self):
        """ plant water potential [cm] """
        return np.array(self.psiXyl)
                    
    def get_es(self, TairC):
        """
            get the atmospheric humidity at saturation (cf. FAO56)
            @param TairC [°C] air temperature
        """
        return 0.61078 * np.exp(17.27 * TairC / (TairC + 237.3))
    
    
    def SPAD2Chl(self,SPAD):
        """ SPAD to leaf chlorophyle content
            ref: Eq 1 Quian 2021, doi: 10.1029/2020JG006076
        @param SPAD [-]
        return leaf chlorophyll content [g cm−2]
        NB: seems to yield too high values (10 times too high?)
        """
        return (0.114 * (SPAD**2) + 7.39 * SPAD + 10.6) * 1e-6
            
    def write_photosynthesis_parameters(self,filename="photosynthesis_parameters"):
        """ write the photosynthesis parameters to a json file """
        m2d_to_cm2s = (24 * 3600) / 1e4 # [m-2 s-1] to [cm-2 d-1]  
        hPa_to_cm = 1.0197
        parameters = {
            "Climate": {
                "Patm": {"value": self.Patm, "unit": "hPa", "description": "Atmospheric pressure"},
                "pCO2": {"value": self.pCO2, "unit": "mol mol-1", "description": "CO2 partial pressure"},
                #"TleafK": {"value": np.array(self.TleafK), "unit": "K", "description": "Leaf temperature"},
                #"Qlight": {"value":np.array( self.Qlight) * (24 * 3600) / 1e4 , "unit": "mol photons cm-2 d-1", "description": "Absorbed photon irradiance"},
                "g_bl": {"value": (np.array( self.g_bl) * m2d_to_cm2s).tolist(), "unit": "mol CO2 cm-2 d-1", "description": "Leaf boundary molar conductance"},
                "g_canopy": {"value": (np.array(self.g_canopy) * m2d_to_cm2s).tolist(), "unit": "mol CO2 cm-2 d-1", "description": "Aerodynamic molar conductance"},
                "g_air": {"value": (np.array(self.g_air) * m2d_to_cm2s).tolist(), "unit": "mol CO2 cm-2 d-1", "description": "Aerodynamic molar conductance"}
            },
            "Plant": {
                "C3 and C4": {
                    "PhotoType": {"value": self.PhotoType,  "description": "Photosynthesis type (C3 or C4)"},
                    "gm": {"value": self.gm * m2d_to_cm2s, "unit": "mol CO2 cm-2 d-1", "description": "Mesophyll resistance"},
                    "g0": {"value": self.g0 * m2d_to_cm2s, "unit": "mol CO2 cm-2 d-1", "description": "Residual stomatal opening"},
                    "a1": {"value": self.a1, "unit": "", "description": "Tuzet 2003 stomatal function parameter"},
                    "alpha": {"value": self.alpha, "unit": "-", "description": "Light response coefficient"},
                    "water stress": {
                        "fwr": {"value": self.fwr, "unit": "-", "description": "Residual opening under water stress"},
                        "fw_cutoff": {"value": self.fw_cutoff, "unit": "-", "description": "To make it easier to get fw"},
                        "sh": {"value": self.sh, "unit": "-", "description": "Sensitivity to water stress"},
                        "p_lcrit": {"value": self.p_lcrit * 1e4 * hPa_to_cm, "unit": "cm", "description": "Min xylem water potential for stomatal opening"}
                    },
                    "Vc": {
                        "Chl": {"value": (np.array(self.Chl) * 1e-6).tolist(), "unit": "g cm-2", "description": "Chlorophyll concentration"},
                        "VcmaxrefChl1": {"value": self.VcmaxrefChl1, "unit": "MPa", "description": "Min psiXil for stomatal opening"},
                        "VcmaxrefChl2": {"value": self.VcmaxrefChl2, "unit": "MPa", "description": "Min psiXil for stomatal opening"}
                    }
                },
                "C3": {
                    "a3": {"value": self.a3, "unit": "-", "description": "maximum electron transport rate to maximum carboxylation rate ratio (Jrefmax = Vcrefmax * a3)"},
                    "theta": {"value": self.theta, "unit": "-", "description": "Curvature of light response"},
                    "oi": {"value": self.oi, "unit": "mol mol-1", "description": "Leaf internal [O2]"}
                },
                "C4": {
                    "Q10_photo": {"value": self.Q10_photo, "unit": "-", "description": "Temperature sensitivity of photosynthesis"},
                }
            }
        }
        
        with open(filename + ".json", "w+") as f:
            json.dump(parameters, f)

    def read_photosynthesis_parameters(self, filename):
        """ read the photosynthesis parameters from a json file """
        m2d_to_cm2s = (24 * 3600) / 1e4 # [m-2 s-1] to [cm-2 d-1] 
        hPa_to_cm = 1.0197
        
        with open(filename + ".json", "r") as f:
            parameters = json.load(f)
            
        # Climate parameters
        self.Patm = parameters["Climate"]["Patm"]["value"]
        self.pCO2 = parameters["Climate"]["pCO2"]["value"]
        self.g_bl = np.array(parameters["Climate"]["g_bl"]["value"]) / m2d_to_cm2s
        self.g_canopy = np.array(parameters["Climate"]["g_canopy"]["value"]) / m2d_to_cm2s
        self.g_air = np.array(parameters["Climate"]["g_air"]["value"]) / m2d_to_cm2s

        # Plant parameters
        self.PhotoType = parameters["Plant"]["C3 and C4"]["PhotoType"]["value"]
        self.gm = parameters["Plant"]["C3 and C4"]["gm"]["value"] / m2d_to_cm2s
        self.g0 = parameters["Plant"]["C3 and C4"]["g0"]["value"] / m2d_to_cm2s
        self.a1 = parameters["Plant"]["C3 and C4"]["a1"]["value"]
        self.alpha = parameters["Plant"]["C3 and C4"]["alpha"]["value"]

        # Water stress parameters
        self.fwr = parameters["Plant"]["C3 and C4"]["water stress"]["fwr"]["value"]
        self.fw_cutoff = parameters["Plant"]["C3 and C4"]["water stress"]["fw_cutoff"]["value"]
        self.sh = parameters["Plant"]["C3 and C4"]["water stress"]["sh"]["value"]
        self.p_lcrit = parameters["Plant"]["C3 and C4"]["water stress"]["p_lcrit"]["value"] / (1e4 * hPa_to_cm)

        # Vc parameters
        self.Chl = np.array(parameters["Plant"]["C3 and C4"]["Vc"]["Chl"]["value"]) / 1e-6
        self.VcmaxrefChl1 = parameters["Plant"]["C3 and C4"]["Vc"]["VcmaxrefChl1"]["value"]
        self.VcmaxrefChl2 = parameters["Plant"]["C3 and C4"]["Vc"]["VcmaxrefChl2"]["value"]

        # C3 parameters
        self.a3 = parameters["Plant"]["C3"]["a3"]["value"]
        self.theta = parameters["Plant"]["C3"]["theta"]["value"]
        self.oi = parameters["Plant"]["C3"]["oi"]["value"]

        # C4 parameters
        self.Q10_photo = parameters["Plant"]["C4"]["Q10_photo"]["value"]
        
        
    def radial_fluxes(self):
        """ plant-exterior exchanges [cm3/day] """
        return np.array(self.outputFlux)