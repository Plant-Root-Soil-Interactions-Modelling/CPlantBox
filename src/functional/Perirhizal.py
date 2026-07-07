import numpy as np
from numpy import linalg as LA
import scipy.linalg as la
from scipy.interpolate import RegularGridInterpolator
from scipy.optimize import fsolve, root_scalar
from scipy.spatial import ConvexHull, Voronoi
from scipy.integrate import odeint
from scipy import integrate
from scipy.special.orthogonal import ts_roots

import time #only used for a test

import math

import plantbox as pb
import plantbox.functional.van_genuchten as vg
from plantbox import MappedSegments, Perirhizal


class PerirhizalPython(Perirhizal):
    """
    Helper class for modelling the perirhizal zone

    Wraps MappedSegments (or specialisations MappedPlant) and adds functions to retrieve information

    * calculates root soil interface potential using steady rate approximation (Schröder et al. 2008)
    * perirhizal_conductance_per_layer calculates the perirhizal conductance in a soil layer (or cell) (Vanderborght et al. 2023, Eqn [6])
    * support of 2D lookup tables (file type is a zipped archive) -> will be IMPROVED
    *
    * analysis across soil grid, e.g. get_density, average, aggregate
    * calculates outer perirhizal radii (based on denisties or Voronoi)

    * TODO maybe seperate geometry related functions (e.g. get_density, get_outer_radii) from soil root interface potential related functions (e.g. soil_root_interface_potentials, create_lookup) into different classes

    run script for different examples, uncomment (to create a 2D lookup table, usage lookup table, and voronoi outer radii)
    """

    def __init__(self, ms=None):
        """@param ms reference to MappedSegments"""
        if ms:
            super().__init__(ms)
        else:
            super().__init__()

        self.lookup_table = None  # optional 2d look up table to find soil root interface potentials
        self.global_lookup_table = None  # optional 3d look up table to find soil root interface potentials
        self.lookup_table_sr_solutes = None # optional 1d look up table for steady state solute flow
        self.lookup_table_sr_solutes_simplified = None # optional 2d lookup tables for the steady rate solute flow (one every diffusion coefficient)
        self.sp = None  # corresponding van genuchten soil parameter
        
        self.alpha_0 = 0.3 # a constant that is used for numerically solving the perirhizal waterflow
        self.h_wilt = -16000
        self.Ds0_ref = 1 # reference diffusion coefficient of solutes in water [cm2/d]
        self.water_filename = "lookup_perirhizal_waterflow_global" #name and location for the global lookup table

    def set_soil(self, sp):
        """sets VG parameters, and no look up table (slow)"""
        vg.create_mfp_lookup(sp)
        self.sp = sp
        self.lookup_table = None
        self.global_lookup_table = None
        self.lookup_table_sr_solutes = None
        self.lookup_table_sr_solutes_simplified = None
 
        
    def open_lookup(self, filename):
        """opens a  look-up table from a file to quickly find soil root interface potentials"""
        npzfile = np.load(filename + ".npz")
        interface = npzfile["interface"]
        inner_kr_b, base_mfp = npzfile["inner_kr_b_"], npzfile["base_mfp_"]
        soil = npzfile["soil"]
        self.lookup_table = RegularGridInterpolator((inner_kr_b, base_mfp), interface)  # method = "nearest" fill_value = None , bounds_error=False
        self.sp = vg.Parameters(soil)
        vg.create_mfp_lookup(self.sp) # I do not know why this has to be repeated here
            
    def open_global_lookup(self, filename):
        """opens a global look-up table from a file to quickly find soil root interface potentials, this lookup table works for all van Genuchten parameter sets"""
        npzfile = np.load(filename + ".npz")
        interface, interface_vg = npzfile["interface"], npzfile["interface_vg"]
        vg_m_, inner_kr_b_, base_mfp_, sx_ = npzfile["vg_m_"], npzfile["inner_kr_b_"], npzfile["base_mfp_"], npzfile["sx_"]
        self.global_lookup_table = RegularGridInterpolator((vg_m_, inner_kr_b_, base_mfp_), interface)  # method = "nearest" fill_value = None , bounds_error=False
        self.lookup_global_mfp = RegularGridInterpolator((vg_m_, sx_), interface_vg)
        #.sp = vg.Parameters(soil)
    
    def open_lookup_solutes_simplified(self, filename):
        """opens a look-up table from a file to quickly find soil root solute concentrations in the steady state case"""
        npzfile = np.load(filename + ".npz")
        integral_AdvDiff_ = npzfile["integral_AdvDiff_"]
        base_mfp_ = npzfile["base_mfp_"]
        soil = npzfile["soil"]
        self.lookup_table_sr_solutes_simplified = RegularGridInterpolator((base_mfp_,[0,1]) , integral_AdvDiff_)  # method = "nearest" fill_value = None , bounds_error=False
        self.sp = vg.Parameters(soil)
        vg.create_mfp_lookup(self.sp) # does this have to be repeated here?
    
    def open_lookup_solutes(self, filename):
        """opens an additional look-up table from a file to quickly find soil root solute concentrations in the steady rate case"""
        npzfile_sr = np.load(filename + ".npz")
        Ds_, outer_mfp_, inner_mfp_ = npzfile_sr["Phi1_"], npzfile_sr["waterflow_"], npzfile_sr["wateruptake_"], npzfile_sr["r_eval_"]
        soil = npzfile_sr["soil"]
        conc_rel_c, inflow_rel_c, Uptake_rel_c = npzfile_sr["conc_rel_c"], npzfile_sr["inflow_rel_c"], npzfile_sr["Uptake_rel_c"]
        conc_mean_c, inflow_mean_c, Uptake_mean_c = npzfile_sr["conc_mean_c"], npzfile_sr["inflow_mean_c"], npzfile_sr["Uptake_mean_c"]
        self.lookup_table_sr_solutes = {
                "conc_rel_c" : [0],
                "conc_mean_c" : [0],
                "inflow_rel_c" : [0],   
                "inflow_mean_c" : [0],
                "Uptake_rel_c" : [0],
                "Uptake_mean_c" : [0]
                }                
        self.lookup_table_sr_solutes["conc_rel_c"] = RegularGridInterpolator((Phi1_, waterflow_, wateruptake_, r_eval_) , conc_rel_c)  # method = "nearest" fill_value = None , bounds_error=False
        self.lookup_table_sr_solutes["conc_mean_c"] = RegularGridInterpolator((Phi1_, waterflow_, wateruptake_, r_eval_) , conc_mean_c)
        self.lookup_table_sr_solutes["inflow_rel_c"] = RegularGridInterpolator((Phi1_, waterflow_, wateruptake_, r_eval_) , inflow_rel_c)
        self.lookup_table_sr_solutes["inflow_mean_c"] = RegularGridInterpolator((Phi1_, waterflow_, wateruptake_, r_eval_) , inflow_mean_c)
        self.lookup_table_sr_solutes["Uptake_rel_c"] = RegularGridInterpolator((Phi1_, waterflow_, wateruptake_, r_eval_) , Uptake_rel_c)
        self.lookup_table_sr_solutes["Uptake_mean_c"] = RegularGridInterpolator((Phi1_, waterflow_, wateruptake_, r_eval_) , Uptake_mean_c)
        self.sp = vg.Parameters(soil)
        vg.create_mfp_lookup(self.sp) # does this have to be repeated here?
        
    def soil_root_interface_potentials(self, rx, sx, inner_kr, rho):
        """
        finds matric potentials at the soil root interface for as all segments
        uses a look up tables if present (see create_lookup, and open_lookup)

        rx             xylem matric potential [cm]
        sx             bulk soil matric potential [cm]
        inner_kr       root radius times hydraulic conductivity [cm/day]
        rho            geometry factor [1] (outer_radius / inner_radius)
        """
        assert len(rx) == len(sx) == len(inner_kr) == len(rho), "rx, sx, inner_kr, and rho must have the same length"
        
        for i in range(len(rx)):
            if rx[i] < -15989:
                rx[i] = -15989
            if sx[i] > -1:
                sx[i] = -1
        
        if self.lookup_table:
            rsx = self.soil_root_interface_potentials_table(rx, sx, inner_kr, rho)
        else:
            if self.global_lookup_table:
                rsx = self.soil_root_interface_potentials_table_global(rx, sx, inner_kr, rho)
            else:
                rsx = np.array([PerirhizalPython.soil_root_interface_(rx[i], sx[i], inner_kr[i], rho[i], self.sp) for i in range(0, len(rx))])
                
        
        return rsx

    @staticmethod
    def soil_root_interface_(rx, sx, inner_kr, rho, sp):
        """
        finds matric potential at the soil root interface for as single segment

        rx             xylem matric potential [cm]
        sx             bulk soil matric potential [cm]
        inner_kr       root radius times hydraulic conductivity [cm/day]
        rho            geometry factor [1] (outer_radius / inner_radius)
        sp             soil parameter: van Genuchten parameter set (type vg.Parameters)
        """
        #vg.create_mfp_lookup(sp)
        #print(vg.fast_mfp[sp](-16000))
        if inner_kr < 1.0e-7:
            return sx
        k_soilfun = lambda hsoil, hint: (vg.fast_mfp[sp](hsoil) - vg.fast_mfp[sp](hint)) * (hsoil - hint) / ((hsoil - hint)**2 + 0.001)  # Vanderborgth et al. 2023, Eqn [7]
        rho2 = np.square(rho)  # rho squared
        b = 2 * (rho2 - 1) / (1 - 0.53 * 0.53 * rho2 + 2 * rho2 * (np.log(rho) + np.log(0.53)))  # Vanderborgth et al. 2023, Eqn [8]
        fun = lambda x: (inner_kr * rx + b * sx * k_soilfun(sx, x)) / (b * k_soilfun(sx, x) + inner_kr) - x
        
        if rx == sx:  # degenerate bracket: no gradient, interface equals bulk potential
            return sx  
        rsx = root_scalar(fun, method="brentq", bracket=[min(rx, sx), max(rx, sx)])
        return rsx.root
     
    @staticmethod
    def soil_root_interface_simp_(self, inner_kr_b, base_mfp, sp):
        """
        finds matric potential at the soil root interface for as single segment
        gives (numerically) the same results as "soil_root_interface_", but is much simpler to create a lookup table for

        rx             xylem matric potential [cm]
        sx             bulk soil matric potential [cm]
        inner_kr       root radius times hydraulic conductivity [cm/day]
        rho            geometry factor [1] (outer_radius / inner_radius)
        sp             soil parameter: van Genuchten parameter set (type vg.Parameters)
        
        inner_kr_b:    inner_kr/b
        base_mfp:      -(inner_kr/b*rx+vg.fast_mfp[sp](sx))
        """
        fun = lambda x: vg.fast_mfp[sp](x) + inner_kr_b * (x - 0*self.h_wilt) +base_mfp
        x_int = [-15999.0,-1.0] # interval for search
        if fun(x_int[0])> 0:
            return x_int[0]
        if fun(x_int[1])< 0:
            return x_int[1]    
        rsx = root_scalar(fun, method="brentq", bracket=[x_int[0], x_int[1]])
        return rsx.root
    
    @staticmethod
    def soil_root_interface_global_(self, inner_kr_b, base_mfp, vg_m):
        """
        finds matric potential at the soil root interface for as single segment

        rx             xylem matric potential [cm]
        sx             bulk soil matric potential [cm]
        inner_kr       root radius times hydraulic conductivity [cm/day]
        rho            geometry factor [1] (outer_radius / inner_radius)
        sp             soil parameter: van Genuchten parameter set (type vg.Parameters)
        
        inner_kr_b:    (inner_kr * alpha)/(alpha_0 * vg_Ks * b)
        base_mfp:      -((inner_kr * alpha)/(alpha_0 * Ks * b) * rx + vg.fast_mfp[sp_base*(alpha_0/alpha)](sx))
        vg_m:          van genuchten parameter m
        
        (sp_base is the van genuchten parameter set of Ks = 1 and alpha = alpha_0, m is an input, theta_r and theta_s are not important for the steady rate model computaiton right now)
        """
        
        
        
        fun = lambda x: self.lookup_global_mfp((vg_m, x)) + inner_kr_b * x +base_mfp
        x_int = [-15990.0,-0.1] # interval for search, get closer to 0 as x = h_sr * alpha / alpha_0, alpha / alpha_0 <= 1
        if fun(x_int[0])> 0:
            return x_int[0]
        if fun(x_int[1])< 0:
            return x_int[1]
        rsx = root_scalar(fun, method="brentq", bracket=[x_int[0], x_int[1]])
        return rsx.root
    
    def solutesuptake_convdiff_(self, watercontent, c_bulk, Vmax, Km, Ds, waterflow, r_root, E, t, sp):
        """
        finds solute concentration at the soil root interface for all segments following T. Roose and G. Kirk 2009 doi:10.1007/s11104-008-9777-z
        
        It assumes a constant watercontent with a constant wateruptake by the root
        
        watercontent   watercontent of the perirhizal zone [cm3/cm3]
        c_bulk         starting solute concentration at the bulk soil [mol/cm3]
        Vmax           maximum solute uptake rate, Michaelis Menten Kinetics [mol/(cm2d)]
        Km             half saturation constant Michaelis Menten Kinetics [mol/cm3]
        Ds             Diffusion constant in water [cm2/d]
        waterflow      steady state waterflow (entire root circumference) [cm2/d]
        r_root         root radius [cm]
        E              minimum net influx into the plant [mol/cm2d]
        t              root age [d]
        sp             van Genuchten parameter set
        """
        assert len(watercontent) == len(c_bulk) == len(Vmax) == len(Km) == len(waterflow) == len(r_root) == len(E) == len(t), "error in Perirhizal.py, solutesuptake_convdiff_: watercontent, c_bulk, Vmax, Km, Ds, waterflow, r_root, E and t must have the same length"
        
        n_segments = len(c_bulk)
        segLength = 1 #reference segment length, should not impact the results [cm]
        
        rsc = np.zeros(n_segments) # solute concentration at the root soil interface
        F = np.zeros(n_segments) # uptake of solutes mol/(cm2d)
        gamma = 0.577 # Eulers constant
        l_func = lambda time : 1/2*np.log(4*np.exp(-gamma)*time+1)
        
        
        for i in range(n_segments):
            #compute diffusion coefficient according to Millington and Quirk
            D = Ds * math.pow(watercontent[i], 10/3) / (sp.theta_S**2)
            
            #unit conversions at the end
            unitconversion_F = Km[i] * D * r_root[i] # to [mol/cm2d] #TODO: look at this again, it doesn't seem right
            unitconversion_c = Km[i] # to [mol/cm3]
            
            #express the inputs without units
            c_inf = c_bulk[i] / Km[i]
            Pe = waterflow[i]/D #Peclet number
            lamb = segLength*r_root[i]*Vmax[i]/(D*Km[i])
            epsilon = segLength*r_root[i]*E[i]/(D*Km[i])
            
            #compute the unitless uptake
            if Pe<=(lamb/(1+c_inf)):
                F[i] = 2*lamb*(c_inf+epsilon*l_func(t[i])) / (1+c_inf+l_func(t[i])*(lamb + epsilon) +math.sqrt(4*(c_inf+epsilon * l_func(t[i]))+(1-c_inf+(lamb-epsilon)*l_func(t[i]))**2)) #Eqn. [9]
            else:
                F[i] = 2*lamb*(c_inf+epsilon*l_func(t[i])) / (1+c_inf+l_func(t[i])*(lamb - Pe + epsilon) +math.sqrt(4*(c_inf+epsilon * l_func(t[i]))*(1-Pe*l_func(t[i]))+(1-c_inf+l_func(t[i])*(lamb-Pe-epsilon))**2)) #Eqm. [10]
            #compute the unitless solute concentration next to the root
            #F(t) = lamb * c / (1+c) - epsilon shortly after Eqn. [10]
            #(1+c)*(F(t)+epsilon) = lamb * c
            #F(t)+epsilon = (-F(t)-epsilon+lamb) * c
            rsc[i] = (F[i]+epsilon)/(lamb-F[i]-epsilon)
            
            #add units
            F[i] = F[i] * unitconversion_F
            rsc[i] = rsc[i] * unitconversion_c
        
        return rsc
    
    def soil_root_solutes_sr(self, Phi_soil, rootwateruptake, waterinflow, r_root, r_prhiz, c_soil, c_outer, Vmax, Km, Ds, sp, mode = "ff"):
        """
        steady rate assumption of solute uptake by roots TODO: insert citaiton
        
        Phi_soil           outer matrix flux potential [cm2/d]
        rootwateruptake    radial root water uptake [cm2/d]
        waterinflow        radial water inflow at r_prhiz [cm2/d]
        r_root             root radius [cm]
        r_prhiz            outer radius [cm]
        c_soil             solute concentration of the cylinder [mol/cm3]
        c_outer            solute concentration outside of the cylinder [mol/cm3], only used in the general steady rate case
        Vmax               maximum solute uptake rate, Michaelis Menten Kinetics [mol/(cm2d)]
        Km                 half saturation constant Michaelis Menten Kinetics [mol/cm3]
        Ds                 Diffusion constant in water [cm2/d]
        sp                 van Genuchten parameter set
        mode               what mode is used for the solute uptake? 
                           "ss" for steady state solute uptake
                           "sr" steady rate solute uptake with a no flux outer BC
                           "ff" for steady rate uptake in combination with the far field approximation to give an outer Cauchy BC #TODO: cite paper for far field approximation
        
        output:
        rsc                solute concentration next to the root [mol/cm3]
        Uptake             solute uptake of the root [mol/(cm d)]
        
        """
        assert len(Phi_soil) == len(rootwateruptake) == len(waterinflow) == len(r_root) == len(r_prhiz) == len(c_soil) == len(c_outer) == len(Vmax) == len(Km), "Phi_soil, rootwateruptake, waterinflow, r_root, r_prhiz, c_soil, c_outer, Vmax and Km must have the same length"
        
        n_segments = len(c_soil)
        
        rsc = np.zeros(n_segments) #conentration next to the root
        Uptake = np.zeros(n_segments) #total solute uptake of the root
        ss_uptake = np.zeros(n_segments) #inflow into the parirhizal zone
        sr_uptake = np.zeros(n_segments) #solute depletion rate of the perirhizal zone
        
        #solve quadratic eqation for root uptake # TODO: Link to publication
        for i in range(n_segments):
            rho = r_prhiz[i] / r_root[i]
            Phi_out =  Phi_soil[i] -( (rootwateruptake[i]-waterinflow[i])*((0.53**2)*(rho**2)/(2*(1-rho**2))+rho**2/(1-rho**2)*(np.log(1/0.53)-0.5)) + waterinflow[i]*np.log(0.53) )
            #Phi = lambda r_rel : Phi_out + (rootwateruptake[i]-waterinflow[i])*((r_rel**2)*(rho**2)/(2*(1-rho**2))+rho**2/(1-rho**2)*(np.log(1/r_rel)-0.5)) + waterinflow[i]*np.log(r_rel) #TODO: test using Phi_soil instead of Phi_outer
            Phi = lambda r : Phi_out + (rootwateruptake[i] - waterinflow[i]) * ((r/r_prhiz[i])**2/2-1/2+np.log(r_prhiz[i]/r)) * rho**2 / (1 - rho**2) + waterinflow[i] * np.log(r / r_prhiz[i])
            Phi = lambda r : 1.0 #TODO: remove times 3
            radial_waterflow = lambda r : (waterinflow[i] + (rootwateruptake[i]-waterinflow[i]) * (r_prhiz[i]**2 - r**2) / (r_prhiz[i]**2 - r_root[i]**2)) / (2 * np.pi * r)
            Ds_func = lambda r_rel : Ds[i] * math.pow(vg.water_content(vg.fast_imfp[sp](Phi(r_rel)),sp),10/3) / (sp.theta_S**2) # Millington and Quirk 
            
            Ds0 = self.Ds0_ref 
            scaling = math.sqrt(Ds[i] / Ds0) #I am an idiot I am an idiot. I messed up the ordering
            Phi1 = Phi(1/scaling) - waterinflow[i] * np.log(scaling) - (rootwateruptake[i] - waterinflow[i]) * r_prhiz[i]**2 / (r_prhiz[i]**2 - r_root[i]**2) * (np.log(scaling)-1/2 * (1-1/scaling**2))
            waterflow = radial_waterflow(scaling) / scaling 
            wateruptake = (rootwateruptake[i] - waterinflow[i]) * r_prhiz[i]**2 / (r_prhiz[i]**2 - r_root[i]**2) / scaling**2 
            
            r_eval = np.logspace(np.log10(r_root[i]), np.log10(r_prhiz[i]), num = 40)
            r_eval = [r * scaling for r in r_eval]
            Phi1 = 1
            waterflow = 0
            wateruptake = 0
            print("r_eval equation before",r_root[i], r_prhiz[i],r_eval)
            conc_rel_c, conc_mean_c, inflow_rel_c, inflow_mean_c, Uptake_rel_c, Uptake_mean_c = self.solute_linearequation_sr_(Phi1, waterflow, wateruptake, r_eval, sp)
            
            print("r_eval equation",r_root[i], r_prhiz[i],r_eval)
            
            #default prefactors for the steady state case
            pre_c = 0
            pre_srUptake = 1
            pre_inflow = 0
            absolute = 0
            #pre_c * c(r_prhiz) + pre_srUptake * srUptake + pre_inflow * inflow = absolute
            
            match mode:
                case "ss": #steady state
                    pre_c, pre_srUptake, pre_inflow, absolute = 0, 1, 0, 0
                case "sr": #steady rate uptake, no influx
                    pre_c, pre_srUptake, pre_inflow, absolute = 0, -1, 1, 0
                case "ff": #far field approximation # TODO insert link to Tiina Roose publication, and current application
                    B_abs = c_outer[i] * (r_prhiz[i]**2)/ (math.exp(-2*r_prhiz[i]**2)-math.exp(-r_prhiz[i]**2))
                    B_cprhiz = -1 * (r_prhiz[i]**2) / (math.exp(-2*r_prhiz[i]**2)-math.exp(-r_prhiz[i]**2))
                    
                    pre_inflow = 1
                    pre_c = - Ds_func(r_prhiz[i]) * B_cprhiz * math.exp(-r_prhiz[i]**2) / r_prhiz[i]**2 - waterinflow[i] / (2 * np.pi * r_prhiz[i]) #TODO check sign
                    absolute = Ds_func(r_prhiz[i]) * B_abs * math.exp(-r_prhiz[i]**2) / r_prhiz[i]**2
                    #pre_c = - Ds_func(r_prhiz[i]) * B_cprhiz * math.exp(-r_prhiz[i]) / r_prhiz[i] - waterinflow[i] / (2 * np.pi * r_prhiz[i]) #TODO check sign
                    #absolute = Ds_func(r_prhiz[i]) * B_abs * math.exp(-r_prhiz[i]) / r_prhiz[i]
                    pre_srUptake = 1
                case "dirichlet":
                    pre_c = 1
                    absolute =  c_outer[i]                   
                case _:    #default
                    pre_c, pre_srUptake, pre_inflow, absolute = 0, 1, 0, 0
            
            #variables to eliminate: ss flow, sr uptake, c_prhiz, Uptake
            # variable for the final equations: rsc
            A = np.zeros((4,4))
            b = np.zeros((4,2)) #right hand side of the linear equation, one for absolute values, one for the rsc value
            
            #linear equations for c_prhiz, rsc, ss_flow, sr_uptake, Uptake (4 linear, 1 quadratic)
            
            #c11 * ss flow + c12 * sr uptake + c13 * rsc = c_mean
            A[0,0] = inflow_mean_c[-1]
            A[0,1] = Uptake_mean_c[-1]
            b[0,1] = -conc_mean_c[-1]
            b[0,0] = c_soil[i]
            
            #c21 * ss flow + c22 * sr uptake + c23 * rsc = c_prhiz
            A[1,0] = inflow_rel_c[-1]
            A[1,1] = Uptake_rel_c[-1]
            A[1,2] = -1
            b[1,1] = -conc_rel_c[-1]
            
            #ss flow + sr uptake = Uptake
            A[2,0] = 1
            A[2,1] = 1
            A[2,3] = -1
            
            #pre_c * c(r_prhiz) + pre_srUptake * srUptake + pre_inflow * inflow (= ssflow) = absolute
            A[3,0] = pre_inflow
            A[3,1] = pre_srUptake
            A[3,2] = pre_c
            b[3,0] = absolute
            
            #m*rsc+n = [ss uptake, sr uptake, c(rprhiz), Uptake]
            n = la.solve(A,b[:,0])
            m = la.solve(A,b[:,1])
            
            #solve quadratic equation 
            #Uptake = Vmax * rsc / (Km + rsc)
            #(Km + rsc) * (m[3]*rsc+n[3]) = Vmax * rsc
            # m[3] * rsc**2 + (Km[i]*m[3]+n[3]-Vmax[i]) * rsc + Km[i]*n[3] = 0 
            
            tol = 1e-9 #arbitrary for now, TODO think of a reasonable tolerance
            
            if abs(m[3])<tol:
                rsc[i] = - (Km[i]*n[3]) / (Km[i]*m[3]+n[3]-Vmax[i])
            else:
                A = m[3] #TODO: use a different letter here than the name of the matrix
                B = (Km[i]*m[3]+n[3]-Vmax[i])
                C = Km[i]*n[3]
                p = B/A
                q = C/A
                if p<0:
                    rsc[i] = - p/2 - math.sqrt(p**2/4-q)
                else:
                    rsc[i] = - p/2 + math.sqrt(p**2/4-q)
            Uptake[i] = m[3]*rsc[i]+n[3]
            ss_uptake[i] = m[0]*rsc[i]+n[0]
            sr_uptake[i] = m[1]*rsc[i]+n[1]
            
            print("output", conc_mean_c[-1], inflow_mean_c[-1], Uptake_mean_c[-1], c_soil)
            print("linear", rsc*conc_mean_c[-1]+ss_uptake* inflow_mean_c[-1]+sr_uptake* Uptake_mean_c[-1], c_soil)
            
        return rsc, Uptake, ss_uptake, sr_uptake, inflow_mean_c
    
    def watersolutes_disc(self, Phi_out, rootwateruptake, waterinflow, r_root, r_prhiz, r_eval, c_root, Ds0, ss_uptake, sr_uptake, sp):
        """
        given the water and solute uptake data this computes the discretisation of the steady rate solutions
        
        Phi_out          outer matrix flux potential [cm2/d]
        rootwateruptake    radial root water uptake [cm2/d]
        waterinflow        radial water inflow at r_prhiz [cm2/d]
        r_root             root radius [cm]
        r_prhiz            outer radius [cm]
        r_eval             positions at which the solute concentration should be evaluated [cm]
        c_root             solute concentration next to the root [mol/cm3]
        Ds0                Diffusion constant in water [cm2/d]
        ss_uptake          outer inflow of solutes [mol/(cm d)]
        sr_uptake          solute uptake from the perirhizal zone [mol/(cm3d)]
        sp                 van Genuchten parameter set
        
        output:
        watercontent, waterpotential, soluteconcentration discretisations
        """
        
        #simple computations
        rho = r_prhiz / r_root
        
        #water
        #Phi_out =  Phi_soil -( (rootwateruptake-waterinflow)*((0.53**2)*(rho**2)/(2*(1-rho**2))+rho**2/(1-rho**2)*(np.log(1/0.53)-0.5)) + waterinflow*np.log(0.53) )
        Phi = lambda r : Phi_out - (rootwateruptake - waterinflow) * ((r/r_prhiz)**2/2-1/2+np.log(r_prhiz/r)) * r_prhiz**2 / (r_prhiz**2 - r_root**2) + waterinflow * np.log(r / r_prhiz)
        Phi = lambda r : 1.0 #TODO: remove times 3
        wateruptake = (rootwateruptake - waterinflow) / (np.pi * (r_prhiz**2 - r_root**2))#TODO doublecheck
        radial_waterflow = lambda r : (waterinflow + wateruptake * np.pi * (r_prhiz**2 - r**2)) / (2 * np.pi * r) #TODO: rename
        waterpotential_func = lambda r : vg.fast_imfp[sp](Phi(r))
        watercontent_func = lambda r : vg.water_content(waterpotential_func(r),sp)
        
        
        #solutes
        #use the subfunction #TODO: alternative lookup table
        scaling = math.sqrt(Ds0 / self.Ds0_ref)
        Phi1 = Phi(1/scaling) - waterinflow * np.log(scaling) - (rootwateruptake - waterinflow) * r_prhiz**2 / (r_prhiz**2 - r_root**2) * (np.log(scaling)-1/2 * (1-1/scaling**2))
        waterflow = radial_waterflow(scaling) / scaling 
        r_eval = [r * scaling for r in r_eval]
        wateruptake = (rootwateruptake - waterinflow) * r_prhiz**2 / (r_prhiz**2 - r_root**2) / scaling**2 
        
        Phi1 = 1
        waterflow = 0
        wateruptake = 0
        
        Vmax = 4.0e-11 * 1 * (2*3.14*r_root) * (24*3600)  #mol/(cm2 s) * cm * cm * (s/d) -> mol / d 
        Vmax_per_area = Vmax / (1 * (2*3.14*r_root)) #mol / d /cm2 = mol/(cm2d)
        Km = 1.5e-7   #mol/cm3
        
        #rsc, Uptake, ss_uptake, sr_uptake, inflow_mean_c_0 = self.soil_root_solutes_sr([Phi_out], [rootwateruptake], [waterinflow], [r_root], [r_prhiz], [1.98e-5], [2.0e-5], [Vmax], [Km], [Ds0], sp, mode = "ss")
        rsc, Uptake, ss_uptake, sr_uptake, inflow_mean_c_0 = self.soil_root_solutes_sr([1], [0], [0], [r_root], [r_prhiz], [1.98e-5], [2.0e-5], [Vmax], [Km], [Ds0], sp, mode = "ss")
        #conc_rel_c, conc_mean_c, inflow_rel_c, inflow_mean_c, Uptake_rel_c, Uptake_mean_c = self.solute_linearequation_sr_(Phi1, waterflow, wateruptake, r_eval, sp)
        conc_rel_c, conc_mean_c, inflow_rel_c, inflow_mean_c, Uptake_rel_c, Uptake_mean_c = self.solute_linearequation_sr_(1, 0, 0, r_eval, sp)
        c_root = rsc
        print("compareinflow", inflow_mean_c_0, inflow_mean_c)
        print("r_eval disc", r_eval)
        
        soluteconcentration = c_root * conc_rel_c + ss_uptake * inflow_rel_c + sr_uptake * Uptake_rel_c   
        soluteconcentration_mean = c_root * conc_mean_c + ss_uptake * inflow_mean_c + sr_uptake * Uptake_mean_c           
        
        print("output disc1", c_root, ss_uptake, sr_uptake, scaling)
        print("output disc2", conc_rel_c[0], inflow_rel_c[0], Uptake_rel_c[0], soluteconcentration[0])
        print("output disc3", conc_rel_c[-1], inflow_rel_c[-1], Uptake_rel_c[-1], soluteconcentration[-1])
        print("output disc4", conc_mean_c[-1], inflow_mean_c[-1], Uptake_mean_c[-1], soluteconcentration_mean[-1])
        
        watercontent = np.array([watercontent_func(r) for r in r_eval])
        waterpotential = np.array([waterpotential_func(r) for r in r_eval])
        
        return watercontent, waterpotential, soluteconcentration
    
    def solute_linearequation_sr_(self, Phi1, waterflow, wateruptake, r_eval, sp):
        """
        computes a linear equation from the diffusion advection ODE on the perirhizal solute flow #TODO: insert equation number
        the radius (and steady rate wateruptake) are scaled by the diffusion coefficient of the solute
        
        Phi1               matrix flux potential at \tilde{r} = 1 [cm2/d]
        waterflow          waterflow at r = 1 [cm/d]
        wateruptake        steady rate wateruptake throughout the perirhizal zone, scaled to \tilde{r} [cm3/(cm3d)]
        r_eval             the equation is solved from r_eval[0] to r_eval[-1], it is a 1d array [cm] 
        
        outputs:
        conc_rel_c         ratio of c(r_rel[i]) to c(r_rel[0]) for d_r c(1)=0
        conc_mean_c        mean solute concentrations for the steady state solute uptake
        inflow_rel_c       ratio of c(r_rel[i]) to d_r c(1) for c(1)=0
        inflow_mean_c      mean solute concentration for c(1)=0
        Uptake_rel_c       ratio of c(r_rel[i]) to rel_soluteUptake for c(1)=0, d_r c(1)=0 (steady rate uptake case)
        Uptake_mean_c      mean solute concentrations relative to  rel_soluteUptake
        """
        
        #reference diffusion coefficient to which r was scaled to [cm2/d]
        Ds0 = self.Ds0_ref 
        
        #r_root = r_eval[0]
        n = len(r_eval)
        
        conc_rel_c = np.ones(n) 
        conc_mean_c = np.ones(n)
        inflow_rel_c = np.ones(n)
        inflow_mean_c = np.ones(n)  
        Uptake_rel_c = np.ones(n)
        Uptake_mean_c = np.ones(n) 

        #used for weighting
        watercontent_times_r = np.ones(n)      
                
        #waterflow
        Phi = lambda r : Phi1 - wateruptake * (r**2/2 - 1/2 - np.log(r)) + waterflow * 2 * np.pi * 1 * np.log(r) # TODO add citation
        Phi = lambda r : 1.0 #TODO: remove times 3
        radial_waterflow = lambda r : (2 * np.pi * 1 * waterflow + wateruptake * np.pi * (1 - r**2)) #TODO: distinct waterflow [cm/d] and radial waterflow [cm2/d]
        watercontent = lambda r : vg.water_content(vg.fast_imfp[sp](Phi(r)), sp)
        Ds = lambda r : Ds0 * math.pow(watercontent(r),10/3) / (sp.theta_S**2)
        
        f_homogen = lambda c, r : ( radial_waterflow(r) * c) / (Ds(r) * r) #TODO: add reference
        f_steadystate = lambda c, r : (1 + radial_waterflow(r) * c) / (Ds(r) * r) #TODO: add reference
        f_quadratic = lambda c, r : (r**2 + radial_waterflow(r) * c) / (Ds(r) * r) #TODO: add reference

        #numerically solve the ODE c' = f(r_rel,c) starting at r_rel=1
        conc_rel_c = odeint(f_homogen, y0 = 1, t = r_eval)
        conc_rel_c = np.array([element[0] for element in conc_rel_c])
        inflow_rel_c = odeint(f_steadystate, y0 = 0, t = r_eval)
        inflow_rel_c = np.array([element[0] for element in inflow_rel_c])
        Uptake_rel_c = odeint(f_quadratic, y0 = 0, t = r_eval)
        Uptake_rel_c = np.array([element[0] for element in Uptake_rel_c])
        
        #compute the means
        conc_mean_c[0] = conc_rel_c[0]
        inflow_mean_c[0] = inflow_rel_c[0]
        Uptake_mean_c[0] = Uptake_rel_c[0]
        watercontent_times_r = np.array([watercontent(r)*r for r in r_eval])
        conc_mean_c[1:] = np.array([np.average(conc_rel_c[0:(j+1)], weights=watercontent_times_r[0:(j+1)]) for j in range(1,n)])      #TODO optiize this computation (the iteration)  
        inflow_mean_c[1:] = np.array([np.average(inflow_rel_c[0:(j+1)], weights=watercontent_times_r[0:(j+1)]) for j in range(1,n)])       
        Uptake_mean_c[1:] = np.array([np.average(Uptake_rel_c[0:(j+1)], weights=watercontent_times_r[0:(j+1)]) for j in range(1,n)])       
        
        return conc_rel_c, conc_mean_c, inflow_rel_c, inflow_mean_c, Uptake_rel_c, Uptake_mean_c
    
    
    def soil_root_solutes_steadyrate_simplified_(self, Phi_root, Phi_soil, r_root, r_prhiz, c_bulk, Vmax, Km, Ds, waterflow, sp, n_approx = 1):
        """
        finds solute concentration at the soil root interface for all segments assuming that the waterflow and soluteflow are proportional throughout the perirhizal zone
        in this simplification the mode of soluteuptake (steady state, steady rate with no influx, other steady rate) is prescribed by the water uptake
        uses a look up tables if present (see create_lookup, and open_lookup) #TODO: add links to the exact lookup table

        Phi_root       matrix flux potential at the root-soil-interface [cm2/d]
        Phi_soil       matrix flux potential at the bulk soil [cm2/d]
        r_root         root radius [cm]
        r_prhiz        perirhizal radius [cm]
        c_bulk         solute concentration at the bulk soil [mol/cm3]
        Vmax           maximum solute uptake rate, Michaelis Menten Kinetics [mol/(cm2d)]
        Km             half saturation constant Michaelis Menten Kinetics [mol/cm3]
        Ds             Diffusion constant in water [cm2/d]
        waterflow      steady state waterflow (entire root circumference) [cm2/d]
        sp             van Genuchten parameter set
        n_approx       how many points are used to approximate the perirhizal solute concentration? 
                       n_approx = 1 means only approximation at the perirhizal radius        
        """
        assert len(Phi_root) == len(Phi_soil) == len(c_bulk) == len(Vmax) == len(Km) == len(waterflow), "Phi_root, Phi_soil, c_bulk, Vmax, Km and waterflow must have the same length"
        
        n_segments = len(c_bulk)
        
        #the radii at which the solute concentration will be tested. They are given as relative values between r_root and r_prhiz. r=1 means r_prhiz.
        # approx(int[fun,0,1]) = sum(weights[i] * fun(x(i)))
        if n_approx ==1:
            x = [1]
            weights = [1]
        else:
            [x,weights] = ts_roots(n_approx) # returns the nodes and weights for the chebyshev numterical integration
            
        
        rho = [r_prhiz[i]/r_root[i] for i in range(n_segments)]
        c_sol_mean2root = np.zeros(n_segments) # ratio of the mean solute concentration to solute concentration next to the root
        mean_watercontent = np.zeros(n_segments) #watercontent next to the root
        rsc = np.zeros(n_segments) #solute concentration next to the root (output)
        F = np.zeros(n_segments) #F is a helper value for computing the ratio between mean solute concentraiton and solute concentration next to the root
        F0 = np.zeros(n_segments) #helper value for the ratio of solute next to the root
        
        #compute a prefactor
        D_tilde = 1/Ds/math.pow(sp.theta_S-sp.theta_R,13/3)*(sp.theta_S*sp.theta_S)
        
        #determine mfp depending on r
        Phi_A = np.zeros(n_segments)
        Phi_C = np.zeros(n_segments)
        for i in range(n_segments):
            Phi_A[i], Phi_C[i] = self.determine_mfp_function(Phi_root[i], Phi_soil[i], rho[i]) #Phi(r/r_prhiz)= A(s^2-ln(s^2))+C


        
        #determine F0
        if self.lookup_table_sr_solutes_simplified:
            F0=[self.lookup_table_sr_solutes_simplified((Phi_root[i],0)) for i in range(n_segments)]
        else:
            F0=[self.integral_AdvectionDiffusion_(Phi_root[i],self.sp) for i in range(n_segments)]
        
        #compute the ratio of (c_mean-C_SW) to (c_root-C_SW) (weighted by water content and cylinder volume)
        for j in range(n_approx):
            test_radii = np.array([math.sqrt(x[j] * (r_prhiz[i]**2-r_root[i]**2)+r_root[i]**2) for i in range(n_segments)]) #scale the chebyshev nodes according to volume
            s = [test_radii[i]/r_prhiz[i] for i in range(n_segments)]
            Phi_current = [Phi_A[i]* (s[i]**2-2*np.log(s[i])) + Phi_C[i] for i in range(n_segments)]
            current_watercontent = [vg.water_content(vg.fast_imfp[self.sp](Phi_current[i]),self.sp) for i in range(n_segments)]
            if self.lookup_table_sr_solutes_simplified:
                F=[(self.lookup_table_sr_solutes_simplified((Phi_current[i],0))-F0[i]) for i in range(n_segments)]
            else:
                F=[(self.integral_AdvectionDiffusion_(Phi_current[i],self.sp)-F0[i]) for i in range(n_segments)]
            for i in range(n_segments):
                c_sol_mean2root[i] += weights[j] * math.exp(D_tilde*F[i]) * current_watercontent[i]
                mean_watercontent[i] += weights[j] * current_watercontent[i]
        for i in range(n_segments):
            c_sol_mean2root[i] = c_sol_mean2root[i] / mean_watercontent[i]
        
        #solve quadratic eqation # TODO: Link to publication
        for i in range(n_segments):
            a1=c_bulk[i]/c_sol_mean2root[i]
            a2=(1-1/c_sol_mean2root[i])/(waterflow[i])
            p=Km[i]-a2*Vmax[i]-a1
            q=-Km[i]*a1
            if p<0:
                rsc[i]=-p/2-math.sqrt(pow(p/2,2)-q)
            else:
                rsc[i]=-p/2+math.sqrt(pow(p/2,2)-q)                
        
        return rsc
    
    @staticmethod
    def integral_AdvectionDiffusion_(Phi_input, sp):
        """
        Integral of (D_s(theta_s-theta_r)^(13/3))/(theta*D(theta)) from Phi = 0.001 to Phi_input

        Phi             matrix flux potential [cm2/d]
        sp              soil parameter: van Genuchten parameter set (type vg.Parameters)
        """
        if Phi_input <=0:
            return 0
        theta_rel = sp.theta_R/(sp.theta_S-sp.theta_R)
        integral_fun = lambda Phi: pow(theta_rel+vg.effective_saturation(vg.fast_imfp[sp](Phi),sp),-13/3)
        integral_AdvDiff, _ = integrate.quad(integral_fun, 1.0e-3, Phi_input)
        
        return integral_AdvDiff
    
        
    def determine_mfp_function(self, Phi_root, Phi_soil, rho):
        """
        determine the spatial function of the mfp
        Phi_root       matrix flux potential at the root-soil-interface [cm2/d]
        Phi_soil       matrix flux potential at the bulk soil [cm2/d]
        rho            geometry factor (outer_radius / inner_radius) [1]
        
        Phi(r/r_prhiz)= A(s^2-ln(s^2))+C
        (steady rate kinetics with a no flux outer BC at s = 1
        
        A (0.53^2-ln(0.53^2))+C=Phi_soil
        A (0.01^2-ln(0.01^2))+C=Phi_root
        
        A = (Phi_soil-Phi_root)/(0.53^2-ln(0.53^2)-0.01^2+ln(0.01^2)
        C = ((0.53^2-ln(0.53^2)) * Phi_root - (0.01^2-ln(0.01^2)) * Phi_soil)/(0.53^2-ln(0.53^2)-0.01^2+ln(0.01^2)
        """ 
        
        det = 0.53**2-2*math.log(0.53)-(1/rho)**2+2*math.log(1/rho)
        a = 0.53**2-2*math.log(0.53)
        c = (1/rho)**2-2*math.log(1/rho)
            
        Phi_A = min((Phi_soil-Phi_root) / det, 0)
        Phi_C = max((a*Phi_root-c*Phi_soil) / det + c * Phi_A, 1.0e-6) - c * Phi_A
        
        return Phi_A, Phi_C


    def perirhizal_conductance_per_layer(self, h_bs, h_sr, sp):
        """
        The perirhizal conductance in a soil layer (or cell) (Vanderborght et al. 2023, Eqn [6]) [day-1],
        see also PlantHydraulicModel.get_soil_rootsystem_concductance().
        If there are no roots in a layer nan is returned

        h_bs           bulk soil matric potential
        h_sr           matric potential at the soil root interface
        sp             soil parameter: van Genuchten parameter set (type vg.Parameters)
        """
        r_root = self.average(self.ms.getEffectiveRadii())  # average radius over layer [cm]
        rld = self.get_density("length")  # root length density per layer [cm-2] = [cm/cm3]
        r_phriz = np.divide(1.0, np.sqrt(np.pi * rld))  # outer perirhizal radii per layer [cm]
        rho = np.divide(r_phriz, r_root)  # [1], see Vanderborght et al. (2023), Eqn [9]
        rho2 = np.square(rho)  # [1]
        b = np.divide(2 * (rho2 - np.ones(rho2.shape)), np.ones(rho2.shape) - 0.53 * 0.53 * rho2 + 2 * rho2 * (np.log(rho) + np.log(0.53)))  # [1], see Eqn [8]
        k_prhiz = lambda h_bs, h_sr: (vg.fast_mfp[sp](h_bs) - vg.fast_mfp[sp](h_sr)) / (h_bs - h_sr)  # an effective conductivity [cm/day]
        dz = (self.ms.maxBound.z - self.ms.minBound.z) / self.ms.resolution.z  # [cm]
        l_root = rld * dz  # [cm-1]
        return 2 * np.pi * np.multiply(l_root, np.multiply(b, k_prhiz(h_bs, h_sr)))  # [day-1], see Vanderborght et al. (2023), Eqn [6]

    # depricated, leave the function in, in case we ever want to compute big lookup tables again
    def create_lookup_mpi(self, filename, sp):
        """
        depricated
        
        Precomputes all soil root interface potentials for a specific soil type
        and saves results into a 4D lookup table

        filename       three files are written (filename, filename_, and filename_soil)
        sp             van genuchten soil parameters, , call
                       vg.create_mfp_lookup(sp) before
        """
        
        if rank == 0:
            print("The function *create_lookup_mpi* is deprecated, use *create_lookup* without mpi instead.")
        
        from mpi4py import MPI

        comm = MPI.COMM_WORLD
        size = comm.Get_size()
        rank = comm.Get_rank()

        rxn = 150
        rx_ = -np.logspace(np.log10(1.0), np.log10(16000), rxn)
        rx_ = rx_ + np.ones((rxn,))
        rx_ = rx_[::-1]
        sxn = 150
        sx_ = -np.logspace(np.log10(1.0), np.log10(16000), sxn)
        sx_ = sx_ + np.ones((sxn,))
        sx_ = sx_[::-1]
        akrn = 100
        akr_ = np.logspace(np.log10(1.0e-7), np.log10(1.0e-4), akrn)
        rhon = 30
        rho_ = np.logspace(np.log10(2.0), np.log10(200.0), rhon)

        if rank == 0:
            print(filename, "calculating", rxn * sxn * rhon * akrn, "supporting points on", size, "thread(s)")

        work_size = rxn * sxn * akrn * rhon
        count = work_size // size  # number of points for each process to analyze
        remainder = work_size % size  # extra points if work_size is not a multiple of size

        if rank < remainder:  # processes with rank < remainder analyze one extra point
            start = rank * (count + 1)  # index of first point to analyze
            stop = start + count + 1  # index of last point to analyze
        else:
            start = rank * count + remainder
            stop = start + count

        interface_local = np.zeros(stop - start)

        # Loop over the range assigned to this rank
        for index_, index in enumerate(range(start, stop)):
            # Convert flat index to multi-dimensional indices
            i = (index // (sxn * akrn * rhon)) % rxn
            j = (index // (akrn * rhon)) % sxn
            k = (index // rhon) % akrn
            l = index % rhon

            if rank == 0 and index_ % 100 == 0:  # to follow progress
                print("at index", index_ + 1, "/", stop - start, "on thread", rank)

            rx = rx_[i]
            sx = sx_[j]
            akr = akr_[k]
            rho = rho_[l]

            interface_local[index_] = PerirhizalPython.soil_root_interface_(rx, sx, akr, rho, sp)

        # data2share needs to use floats for Allgatherv with MPI.DOUBLE to work.
        data2share = np.array(interface_local, dtype=np.float64)

        # other data needed by comm.Allgatherv
        all_sizes = np.full(size, count)
        all_sizes[:remainder] += 1

        offsets = np.zeros(len(all_sizes), dtype=np.int64)
        offsets[1:] = np.cumsum(all_sizes)[:-1]
        all_sizes = tuple(all_sizes)
        offsets = tuple(offsets)

        # share the vectors
        interface = np.zeros(work_size)
        comm.Allgatherv([data2share, MPI.DOUBLE], [interface, all_sizes, offsets, MPI.DOUBLE])

        if rank == 0:
            interface = interface.reshape(rxn, sxn, akrn, rhon)  # reset shape
            np.savez(filename, interface=interface, rx_=rx_, sx_=sx_, akr_=akr_, rho_=rho_, soil=list(sp))
            self.lookup_table = RegularGridInterpolator((rx_, sx_, akr_, rho_), interface)
            self.sp = sp

    def create_lookup(self, filename, sp):
        """
        Precomputes all soil root interface potentials for a specific soil type
        and saves results into a 2D lookup table

        filename       three files are written (filename, filename_, and filename_soil)
        sp             van genuchten soil parameters, , call
                       vg.create_mfp_lookup(sp) before
        """
        
        
        # intervals for all inputs
        rx_int_abs = [1.0, 15999.0] #absolute values for logarithmic scaling
        sx_int_abs = [1.0, 15999.0]
        akr_int = [1.0e-7,1.0e-4]
        rho_int = [2.0,200.0]
        
        #the lookup table only needs 2 inputs from combinations of the 4 general inputs
        b_func = lambda rho: 2 * (rho*rho - 1) / (1 - 0.53 * 0.53 * rho*rho + 2 * rho*rho * (np.log(rho) + np.log(0.53)))  # Vanderborgth et al. 2023, Eqn [8]
        b_0 = b_func(rho_int[0])
        b_1 = b_func(rho_int[1])
        b_int = [min(b_0,b_1),max(b_0,b_1)]
        
        #compute the intervals of these 2 possible inputs for the lookup table
        
        inner_kr_bn = 100
        inner_kr_b_int = [akr_int[0]/b_int[1],akr_int[1]/b_int[0]]
        inner_kr_b_= np.logspace(np.log10(inner_kr_b_int[0]), np.log10(inner_kr_b_int[1]), inner_kr_bn)
        
        min_int, max_int = vg.fast_mfp[sp](-sx_int_abs[1]), vg.fast_mfp[sp](-sx_int_abs[0])
        base_mfp_n_1 = 500
        base_mfp_n_2 = 500
        base_mfp_n = base_mfp_n_1 + base_mfp_n_2
        base_mfp_int = [inner_kr_b_int[0]*(0*self.h_wilt-(-rx_int_abs[0])) - vg.fast_mfp[sp](-sx_int_abs[0]), inner_kr_b_int[1]*(0*self.h_wilt-(-rx_int_abs[1])) - vg.fast_mfp[sp](-sx_int_abs[1])]
        #base_mfp_int are negative and positive
        base_mfp_interval1= -np.logspace(np.log10(-base_mfp_int[0]), np.log10(1.0e-9), base_mfp_n_1)
        base_mfp_interval2= np.logspace(np.log10(1.0e-9), np.log10(base_mfp_int[1]), base_mfp_n_2)
        base_mfp_ = np.concatenate((base_mfp_interval1,base_mfp_interval2))
        
        interface = np.zeros((inner_kr_bn,base_mfp_n))
        
        print("Creating a lookup table for the hydraulic perirhizal resistance model")
        
        for i, inner_kr_b in enumerate(inner_kr_b_):
            print(i, " / ", inner_kr_bn)
            for j, base_mfp in enumerate(base_mfp_):
                interface[i, j] = PerirhizalPython.soil_root_interface_simp_(self, inner_kr_b, base_mfp, sp)
        np.savez(filename, interface=interface, inner_kr_b_=inner_kr_b_, base_mfp_=base_mfp_, soil=list(sp))
        self.lookup_table = RegularGridInterpolator((inner_kr_b_, base_mfp_), interface)
        self.sp = sp
        
        print("Done with creating a lookup table for the hydraulic perirhizal resistance model")
        return 0

    def create_lookup_global(self, filename, dummy_sp = vg.Parameters([0.078, 0.43, 0.036, 1.56, 24.96])):
        """
        Precomputes all soil root interface potentials for a specific soil type
        and saves results into a 3D lookup table

        filename       three files are written (filename, filename_, and filename_soil)
        dummy_sp       van genuchten soil parameters set for the perirhizal zone, not needed in the actual lookup table
                       has been set to hydrus loam as a default
        """

        
        # intervals for all inputs
        rx_int_abs = [0.1, 16000.0] #absolute values for logarithmic scaling
        sx_int_abs = [0.09, 16000.0]
        akr_int = [1.0e-7,1.0e-4]
        rho_int = [2.0,200.0]
        vg_m_int = [0.1,1.5]
        vg_Ks_int = [1.0,100.0]
        vg_alpha_int = [0.01,0.3]
        
        #the lookup table only needs 2 inputs from combinations of the 4 general inputs
        b_func = lambda rho: 2 * (rho*rho - 1) / (1 - 0.53 * 0.53 * rho*rho + 2 * rho*rho * (np.log(rho) + np.log(0.53)))  # Vanderborgth et al. 2023, Eqn [8]
        b_0 = b_func(rho_int[0])
        b_1 = b_func(rho_int[1])
        b_int = [min(b_0,b_1),max(b_0,b_1)]
        
        vg_m_n = 100 
        vg_m_ = np.logspace(np.log10(vg_m_int[0]), np.log10(vg_m_int[1]), vg_m_n)
        
        sx_n = 200
        sx_ = -np.logspace(np.log10(sx_int_abs[0]), np.log10(sx_int_abs[1]), sx_n) 
        
        interface_vg = np.zeros((vg_m_n,sx_n))
        
        #construct a smaller lookup table for the simplified matrix flux potential
        print("Creating a general lookup table for the matrix flux potential")
        for i, vg_m in enumerate(vg_m_):
            print(i, " / ", vg_m_n)
            sp_dummy = vg.Parameters([0.1,0.4,self.alpha_0,1./(1.-vg_m),1.0])
            vg.create_mfp_lookup(sp_dummy, verbose=False)
            for j, sx in enumerate(sx_):
                interface_vg[i, j] = vg.matric_flux_potential(sx, sp_dummy)
        self.lookup_global_mfp = RegularGridInterpolator((vg_m_, sx_), interface_vg)        


        inner_kr_bn = 200
        inner_kr_b_abs_int = [akr_int[0]/(b_int[1]*vg_Ks_int[1]),akr_int[1]/(b_int[0]*vg_Ks_int[0])]
        inner_kr_b_= np.logspace(np.log10(inner_kr_b_abs_int[0]), np.log10(inner_kr_b_abs_int[1]), inner_kr_bn)
        
        #base_mfp_n = 200
        base_mfp_n_1 = 100
        base_mfp_n_2 = 100
        base_mfp_n = base_mfp_n_1 + base_mfp_n_2
        #base_mfp_int are negative and positive
        base_mfp_int = [inner_kr_b_abs_int[0]*(-(-rx_int_abs[0]) * vg_alpha_int[0]/self.alpha_0) - self.lookup_global_mfp((vg_m_int[1],min(-sx_int_abs[0] * vg_alpha_int[0]/self.alpha_0,-1))), inner_kr_b_abs_int[1]*(-(-rx_int_abs[1]) *vg_alpha_int[1]/self.alpha_0) -self.lookup_global_mfp((vg_m_int[0],-sx_int_abs[1] * vg_alpha_int[1]/self.alpha_0))] 
        base_mfp_interval1= -np.logspace(np.log10(-base_mfp_int[0]), np.log10(1.0e-9), base_mfp_n_1)
        base_mfp_interval2= np.logspace(np.log10(1.0e-9), np.log10(base_mfp_int[1]), base_mfp_n_2)
        base_mfp_ = np.concatenate((base_mfp_interval1,base_mfp_interval2))
        
        interface = np.zeros((vg_m_n,inner_kr_bn,base_mfp_n))
        
        tol = 0.001 #tolerance for extreme values of the matrix flux potential at the soil -> output h_rs = h_soil or h_x
        
        print("Creating a global lookup table for the hydraulic perirhizal resistance model")
        for i, vg_m in enumerate(vg_m_):
            print(i, " / ", vg_m_n)
            for j, inner_kr_b in enumerate(inner_kr_b_):
                for k, base_mfp in enumerate(base_mfp_):
                    interface[i, j, k] = PerirhizalPython.soil_root_interface_global_(self, inner_kr_b, base_mfp, vg_m)
        np.savez(filename, interface=interface, interface_vg=interface_vg, vg_m_=vg_m_, inner_kr_b_=inner_kr_b_, base_mfp_=base_mfp_, sx_=sx_, soil=list(sp))
        self.global_lookup_table = RegularGridInterpolator((vg_m_, inner_kr_b_, base_mfp_), interface)
        self.sp = dummy_sp
        
        print("Done with creating a general lookup table for the matrix flux potential")
        return 0
    
    def create_integralAdvectionDiffusion_lookup(self, filename, sp):
        """      
        Precomputes all integrals for the steady state solute flow
        
        Phi       upper matrix flux potential (bottom is set to 0) [cm2/d]
        sp             van genuchten soil parameters, , call 
                       vg.create_mfp_lookup(sp) before 
        """
        Phin = 300
        Phi_ = np.logspace(np.log10(1.0e-9), np.log10(500), Phin)
        base_mfp_=Phi_
        integral_overD_ = np.zeros((Phin, 2))
        print("Creating a lookup table for the steady state solute flow")
        for i, Phi in enumerate(Phi_):
            print(i, " / ", Phin)
            integral_overD_[i,0] = PerirhizalPython.integral_overDiffusion_(Phi, sp)
            integral_overD_[i,1] = integral_overD_[i,0]
        np.savez(filename, integral_AdvDiff_ = integral_AdvDiff_, base_mfp_ = base_mfp_, soil = list(sp))
        self.lookup_table_solutes = RegularGridInterpolator((base_mfp_,[0,1]) , integral_overD_)
        self.sp = sp
        return integral_overD_, base_mfp_
    
    def create_integralconcentration_lookup(self, filename, Ds_, sp):
        """      
        Precomputes all integrals for the steady state solute flow
        
        Phi       upper matrix flux potential (bottom is set to 0) [cm2/d]
        Ds_        array of all diffusion coefficients [cm2/d]
        sp             van genuchten soil parameters, , call 
                       vg.create_mfp_lookup(sp) before 
        """
        
        #repeat steady state case
        integral_overD_, base_mfp_ = self.create_integralDiffusion_lookup(filename, self.sp)
        
        n_Ds = len(Ds_)
        
        Phin = 200
        Phi_ = np.logspace(np.log10(1.0e-6), np.log10(300), Phin)
        outer_mfp_=Phi_
        inner_mfp_=Phi_
        R_sr_ = np.zeros((n_Ds,Phin,Phin))
        print("Creating a lookup table for the steady rate solute flow")
        for i, Ds in enumerate(Ds_):
            print("starting with diffusion coefficient number ", str(i+1), ":")
            C_d = 1/Ds/math.pow(sp.theta_S-sp.theta_R,13/3)*(sp.theta_S*sp.theta_S)
            for j, Phi_in in enumerate(inner_mfp_):
                print(j)
                for k, Phi_out in enumerate(outer_mfp_):
                    R_sr_[i,j,k] = self.integral_overconcentration_(C_d, Phi_in, Phi_out, 100.0, sp)
        
        np.savez(filename, integral_overD_=integral_overD_, base_mfp_=base_mfp_, Ds_ = Ds_, outer_mfp_ = outer_mfp_, inner_mfp_ = inner_mfp_, R_sr_ = R_sr_, soil = list(sp))
        self.lookup_table_sr_solutes = RegularGridInterpolator((Ds_, inner_mfp_, outer_mfp_) , R_sr_)
        self.sp = sp
    
    def soil_root_interface_potentials_table(self, rx, sx, inner_kr_, rho_):
        """
        finds potential at the soil root interface using a lookup table

        rx             xylem matric potential [cm]
        sx             bulk soil matric potential [cm]
        inner_kr       root radius times hydraulic conductivity [cm/day]
        rho            geometry factor [1]
        f              function to look up the potentials
        """
        
        """
        the original equations were
        vg.fast_mfp[sp](x) = int_(h_wilting)^(x)K(S(h))dh    
        k_soilfun(sx, h_sr)= (vg.fast_mfp[sp](h_sr)-vg.fast_mfp[sp](sx))/(h_sr-sx)  Vanderborgth et al. 2023, Eqn [7]
        (inner_kr * rx + b * sx * k_soilfun(sx, x)) / (b * k_soilfun(sx, x) + inner_kr) - x = 0
         
        equivalent:
        vg.fast_mfp[sp](x) + (inner_kr / b) * x - (inner_kr / b * rx + vg.fast_mfp[sp](sx)) = 0
        vg.fast_mfp[sp](x) + inner_kr_b     * x + base_mfp = 0
        x only depends on 2 variables

        """
        
        try:
            sx = np.array(sx)
            mask = inner_kr_ == 0
            inner_kr_[mask] = 1.0e-7
            
            #compute two inputs for the lookup table
            rho = np.array(rho_)
            rho2 = np.multiply(rho,rho)
            b = 2 * (rho2 - 1) / (1 - 0.53 * 0.53 * rho2 + 2 * rho2 * (np.log(rho) + np.log(0.53)))  # Vanderborgth et al. 2023, Eqn [8]
            inner_kr_b = np.divide(inner_kr_,b)
            base_mfp = - (inner_kr_b * rx + vg.fast_mfp[self.sp](sx))
            
            rsx = self.lookup_table((inner_kr_b, np.array([max(min(base_mfp[i],500),-500) for i in range(len(base_mfp))])))
            rsx[mask] = sx[mask]  # if inner_kr is zero, there is no flow, and the interface potential is the same as the soil potential
        except:
            if np.max(rx) > 0:
                print("xylem matric potential positive", np.max(rx), "at", np.argmax(rx))
            if np.min(rx) < -16000:
                print("xylem matric potential under -16000 cm", np.min(rx), "at", np.argmin(rx))
            if np.max(sx) > 0:
                print("soil matric potential positive", np.max(sx), "at", np.argmax(sx))
            if np.min(sx) < -16000:
                print("soil matric potential under -16000 cm", np.min(sx), "at", np.argmin(sx))
            if np.min(inner_kr_) < 1.0e-7:
                print("radius times radial conductivity below 1.e-7", np.min(inner_kr_), "at", np.argmin(inner_kr_))
            if np.max(inner_kr_) > 1.0e-4:
                print("radius times radial conductivity above 1.e-4", np.max(inner_kr_), "at", np.argmax(inner_kr_))
            if np.min(rho_) < 1:
                print("geometry factor below 1", np.min(rho_), "at", np.argmin(rho_))
            if np.max(rho_) > 200:
                print("geometry factor above 200", np.max(rho_), "at", np.argmax(rho_))

            print("PerirhizalPython.soil_root_interface_potentials_table(): table look up failed, value exceeds table")
            print(
                "\trx",
                np.min(rx),
                np.max(rx),
                "sx",
                np.min(sx),
                np.max(sx),
                "inner_kr",
                np.min(inner_kr_),
                np.max(inner_kr_),
                "rho",
                np.min(rho_),
                np.max(rho_),
            )
            raise

        return rsx
        
    def soil_root_interface_potentials_table_global(self, rx, sx, inner_kr_, rho_):
        """
        finds potential at the soil root interface using a globla lookup table, the lookup table works for all van Genuchten parameter sets

        rx             xylem matric potential [cm]
        sx             bulk soil matric potential [cm]
        inner_kr       root radius times hydraulic conductivity [cm/day]
        rho            geometry factor [1]
        f              function to look up the potentials
        """
        
        """
        the original equations were
        vg.fast_mfp[sp](x) = int_(h_wilting)^(x)K(S(h))dh    
        k_soilfun(sx, h_sr)= (vg.fast_mfp[sp](h_sr)-vg.fast_mfp[sp](sx))/(h_sr-sx)  Vanderborgth et al. 2023, Eqn [7]
        (inner_kr * rx + b * sx * k_soilfun(sx, x)) / (b * k_soilfun(sx, x) + inner_kr) - x = 0
        
        equivalent:
        vg.fast_mfp[sp](x) = alpha_0/alpha*vg_Ks * int_(h_wilting*(alpha/alpha_0))^(x*(alpha/alpha_0))K(S(alpha_0/alpha*h))/vg_Ks dh
        vg.fast_mfp[sp](x) = alpha_0/alpha*vg_Ks * (vg.fast_mfp[sp_base](x*(alpha/alpha_0))-vg.fast_mfp[sp_base](h_wilt*(alpha/alpha_0)))
        vg.fast_mfp[sp](x)-vg.fast_mfp[sp](sx) = alpha_0/alpha*vg.Ks * (vg.fast_mfp[sp_base](x*(alpha/alpha_0))-vg.fast_mfp[sp_base](sx*(alpha/alpha_0)))
        sp_base is the van genuchten parameter set of Ks = 1 and alpha = alpha_0. Since sx*(alpha/alpha_0) should ideally be above the wilting point, choose alpha_0 large
        
        vg.fast_mfp[sp_base](x) + (inner_kr * alpha)/(alpha_0 * vg_Ks * b) * x - ((inner_kr * alpha)/(alpha_0 * Ks * b) * rx + vg.fast_mfp[sp_base](sx)) = 0
        vg.fast_mfp[sp_base](x) + inner_kr_b                               * x + base_mfp = 0

        x only depends on 3 variables: inner_kr_b, base_mfp, vg_m (vg.fast_mfp[sp_base](x) depends only on x and this van Genuchten parameter)
        
        note: this function might not be particular fast as the creation of the sp_base parameter set might take a while, 

        """
        
        try:
            sx = np.array(sx)
            mask = inner_kr_ == 0
            inner_kr_[mask] = 1.0e-7
            
            
            #compute two inputs for the lookup table
            rho = np.array(rho_)
            rho2 = np.multiply(rho,rho)
            b = 2 * (rho2 - 1) / (1 - 0.53 * 0.53 * rho2 + 2 * rho2 * (np.log(rho) + np.log(0.53)))  # Vanderborgth et al. 2023, Eqn [8]
            inner_kr_b = np.divide(inner_kr_,self.sp.Ksat*b)
            base_mfp = - ( inner_kr_b * rx * self.sp.alpha / self.alpha_0 + self.lookup_global_mfp((self.sp.m,sx*self.sp.alpha/self.alpha_0)))
            #lookup table
            rsx = np.array([self.global_lookup_table((self.sp.m, inner_kr_b[i], base_mfp[i])) for i in range(len(inner_kr_b))])*self.alpha_0/self.sp.alpha
            rsx[mask] = sx[mask]  # if inner_kr is zero, there is no flow, and the interface potential is the same as the soil potential
        except:
            if np.max(rx) > 0:
                print("xylem matric potential positive", np.max(rx), "at", np.argmax(rx))
            if np.min(rx) < -16000:
                print("xylem matric potential under -16000 cm", np.min(rx), "at", np.argmin(rx))
            if np.max(sx) > 0:
                print("soil matric potential positive", np.max(sx), "at", np.argmax(sx))
            if np.min(sx) < -16000:
                print("soil matric potential under -16000 cm", np.min(sx), "at", np.argmin(sx))
            if np.min(inner_kr_) < 1.0e-7:
                print("radius times radial conductivity below 1.e-7", np.min(inner_kr_), "at", np.argmin(inner_kr_))
            if np.max(inner_kr_) > 1.0e-4:
                print("radius times radial conductivity above 1.e-4", np.max(inner_kr_), "at", np.argmax(inner_kr_))
            if np.min(rho_) < 1:
                print("geometry factor below 1", np.min(rho_), "at", np.argmin(rho_))
            if np.max(rho_) > 200:
                print("geometry factor above 200", np.max(rho_), "at", np.argmax(rho_))

            print("PerirhizalPython.soil_root_interface_potentials_table(): table look up failed, value exceeds table")
            print(
                "\trx",
                np.min(rx),
                np.max(rx),
                "sx",
                np.min(sx),
                np.max(sx),
                "inner_kr",
                np.min(inner_kr_),
                np.max(inner_kr_),
                "rho",
                np.min(rho_),
                np.max(rho_),
            )
            raise

        return rsx    

    def get_density(self, type: str, volumes=None):
        """retrieves length, surface or volume density [cm/cm3, cm2/cm3, or cm3/cm3] per soil cells

        type       length, surface or volume
        volumes    soil volume cells, if empty a recangular grid is assumed (its data stored in MappedSegments), see get_default_volumes_()
        """
        if type == "length":
            f = lambda a, l: l  # cm
        elif type == "surface":
            f = lambda a, l: l * 2 * a * np.pi  # cm2
        elif type == "volume":
            f = lambda a, l: l * a * a * np.pi  # cm3
        else:
            print("PerirhizalPython.get_density() unknown type (should be 'length', 'surface', or 'volume')", type)
            raise
        if volumes is None:
            volumes = self.get_default_volumes_()
        ms = self.ms  # rename
        cell2seg = ms.cell2seg
        l = ms.segLength()
        radii = ms.getEffectiveRadii()
        n = int(ms.resolution.x * ms.resolution.y * ms.resolution.z)
        d = np.zeros((n,))
        for i in range(0, n):
            if i in cell2seg:
                for si in cell2seg[i]:
                    d[i] += f(radii[si], l[si])
        return np.divide(d, volumes)  # cm/cm3, cm2/cm3, or cm3/cm3

    def aggregate(self, seg_param):
        """sums up params over the grid cells (parrams is given per segment)"""
        ms = self.ms  # rename
        cell2seg = ms.cell2seg  # rename
        n = int(ms.resolution.x * ms.resolution.y * ms.resolution.z)
        d = np.zeros((n,))
        for i in range(0, n):
            if i in cell2seg:
                for si in cell2seg[i]:
                    d[i] += seg_param[si]
        return d

    def average(self, seg_param):
        """averages params over the grid cells (parrams is given per segment)"""
        ms = self.ms  # rename
        cell2seg = ms.cell2seg  # rename
        n = int(ms.resolution.x * ms.resolution.y * ms.resolution.z)
        d = np.zeros((n,))
        for i in range(0, n):
            if i in cell2seg:
                c = len(cell2seg[i])
                for si in cell2seg[i]:
                    d[i] += seg_param[si] / c
        return d

    def to_range_(self, x: list, other=None, min_: float = -np.inf, max_: float = np.inf):
        """returns ths list with x values within min_ and max_ (and drops nans)
        drops the same entries from @param other list
        TODO revise
        """
        y, z = [], []
        for i, x_ in enumerate(x):
            if (not np.isnan(x_)) and x_ >= min_ and x_ <= max_:
                y.append(x_)
                if other is not None:
                    z.append(other[i])
        return y, z

    def get_outer_radii(self, type: str, volumes: list = None):
        """retrieves the outer radii of the perirhizal zones [cm]

        type       each soil volume is splitted into perirhizal zones proportional to segment "length", "surface" or "volume"
        volumes    soil volume cells, if empty a recangular grid is assumed (its data stored in MappedSegments), see get_default_volumes_()
        """
        if type == "length":
            return self.get_outer_radii_(type, volumes)  # Python
            # return super().segOuterRadii(2, volumes) # C++
        elif type == "surface":
            return self.get_outer_radii_(type, volumes)  # Python
            # return super().segOuterRadii(1, volumes) # C++
        elif type == "volume":
            return self.get_outer_radii_(type, volumes)  # Python
            # return super().segOuterRadii(0, volumes) # C++
        elif type == "voronoi" or type == "voronoi_bounded":
            if volumes is not None:
                print("PerirhizalPython.get_outer_radii() Warning: parmeter 'volumes' is set but not applicable for voronoi")
            return self.get_outer_radii_voronoi("bounded")  # default is bounded voronoi
        elif type == "voronoi_periodic":
            if volumes is not None:
                print("PerirhizalPython.get_outer_radii() Warning: parmeter 'volumes' is set but not applicable for voronoi")
            return self.get_outer_radii_voronoi("periodic")
        else:
            print("PerirhizalPython.get_outer_radii() unknown type (should be 'length', 'surface', or 'volume')", type)
            raise

    def get_outer_radii_(self, type: str, volumes: list = None):
        """python version of segOuterRadii, see get_outer_radii()"""
        if type == "length":
            f = lambda a, l: l  # cm
        elif type == "surface":
            f = lambda a, l: l * 2 * a * np.pi  # cm2
        elif type == "volume":
            f = lambda a, l: l * a * a * np.pi  # cm3
        else:
            print("PerirhizalPython.get_outer_radii_() unknown type (should be 'length', 'surface', or 'volume')", type)
            raise
        ms = self.ms  # rename that
        cell2seg = ms.cell2seg
        radii = ms.getEffectiveRadii()
        lengths = ms.segLength()
        width = ms.getDomainWidth()
        outer_r = np.zeros((len(radii),))
        if volumes is None:
            n = ms.resolution.x * ms.resolution.y * ms.resolution.z
            volumes = np.ones(int(n)) * (width.x * width.y * width.z) / n
            # print("PerirhizalPython.get_outer_radii_(): each soil volume has", (width.x * width.y * width.z) / n, "cm3")
        # print("volumes", volumes)
        # print("cell2seg.keys()", cell2seg.keys())
        # print("len(cell2seg[-1])", len(cell2seg[-1]))
        # print("len(cell2seg[0])", len(cell2seg[0]))
        for cell_id, seg_ids in cell2seg.items():
            if cell_id >= 0:
                tt = np.sum(np.array([f(radii[i], lengths[i]) for i in seg_ids]))
                for si in seg_ids:
                    t = f(radii[si], lengths[si]) / tt  # proportionality factor (must sum up to == 1 over cell)
                    v = t * volumes[cell_id]  # target volume
                    outer_r[si] = np.sqrt(v / (np.pi * lengths[si]) + radii[si] * radii[si])

        return outer_r

    def get_outer_radii_voronoi(self, type="bounded"):
        """retrives the outer radii  [cm] and corresponding perirhizal volumes [cm3]
        based on a 3D Voronoi diagram that are bounded by the soil volumes
        using the approach of mirrored nodes (for bounded) or shifted (for periodic).

        TODO doc
        """
        ms = self.ms  # rename
        ns = ms.getNumberOfMappedSegments()
        nodes_ = ms.nodes
        nodes = np.array([[n.x, n.y, n.z] for n in nodes_])  # to numpy array
        # nodes = nodes[: ns + 1]  # TODO there are more nodes than segments,
        # e.g. some nodes are not connected to the root system.For now we just drop them, but we should check why they are there in the first place.
        print("nodes", len(nodes_))
        print("nodes truncated", nodes.shape)
        print("segments", ns)

        x = int(ms.resolution.x)
        y = int(ms.resolution.y)
        z = int(ms.resolution.z)

        width = ms.getDomainWidth()
        # print("making periodic", width)
        # print("nodes", nodes.shape)
        # print("resolution", x, y, z)
        nodes = self.make_periodic_(nodes, width.x, width.y)

        vol = np.zeros((nodes.shape[0]))
        for i in range(0, x):
            for j in range(0, y):
                for k in range(0, z):

                    # print("\nIteration", i, j, k)
                    if type == "bounded":
                        n, ni = self.mirror_(nodes, i, j, k)
                    elif type == "periodic":
                        n, ni = self.shift_(nodes, i, j, k)
                    else:
                        print("PerirhizalPython.get_outer_radii_() unknown type (should be 'bounded', or 'periodic')", type)
                        raise
                    nn = ni.shape[0]
                    print("number of nodes in cell", nn)

                    if nn > 0:
                        vor = Voronoi(n)

                        # print("Number of regions", len(vor.point_region))
                        for idx, reg_num in enumerate(vor.point_region):
                            indices = vor.regions[reg_num]
                            i_ = reg_num - 1
                            if idx < nn:
                                if -1 in indices:  # some regions can be opened
                                    # print(i_)
                                    # print(ni)
                                    vol[ni[idx]] = np.nan
                                    print("nan encountered")
                                else:
                                    vol[ni[idx]] = ConvexHull(vor.vertices[indices]).volume
                                    # print(ni[i_])
                                    # print("index", ni[i_], "vol", vol[ni[i_]])
                            # else:
                            #     if i_ < 0:
                            #         print("When does that happen, again?", idx, nn, indices)

        if np.sum(np.isnan(vol)) > 1:
            print("nan encountered", np.sum(np.isnan(vol)), np.where(np.isnan(vol)))

        outer_r = np.zeros((vol.shape[0] - 1,))
        radii = self.ms.getEffectiveRadii()
        lengths = self.ms.segLength()
        for i, v in enumerate(vol[1:]):  # seg_index = node_index -1
            if v > 0:
                outer_r[i] = np.sqrt(v / (np.pi * lengths[i]) + radii[i] * radii[i])

        return outer_r

    def get_node_mesh(self, nodes):
        """Creates a VTK unstructured grid containing nodes"""
        import vtk

        grid = vtk.vtkUnstructuredGrid()
        points_array = vtk.vtkPoints()
        for n in nodes:
            points_array.InsertNextPoint(n)
        grid.SetPoints(points_array)
        return grid

    def get_voronoi_mesh(self, crop_domain=None):
        """Creates a VTK grid containing vor_nodes; assumes periodic boundaries in x and y direction, and open in z direction"""
        nodes_ = self.ms.nodes  # make nodes periodic and pass all
        nodes = np.array([[n.x, n.y, n.z] for n in nodes_])  # to numpy array
        min_b = self.ms.minBound
        max_b = self.ms.maxBound
        width = np.array([max_b.x, max_b.y, max_b.z]) - np.array([min_b.x, min_b.y, min_b.z])
        nodes = self.make_periodic_(nodes, width[0], width[1])
        return self.get_voronoi_mesh_(nodes, len(nodes), crop_domain)

    def get_voronoi_mesh_(self, vor_nodes, noc, crop_domain=None):
        """Creates a VTK grid containing vor_nodes;
        only cells inside the domain crop_domain are considered (default is bounding box of mapped segments)
        TODO revise
        """
        import vtk

        print("\nget_voronoi_mesh_")

        if not crop_domain:  # always crop with min_b, max_b
            min_b = self.ms.minBound
            max_b = self.ms.maxBound
            crop_domain = pb.SDF_Cuboid(min_b, max_b)

        vor = Voronoi(vor_nodes)
        grid = vtk.vtkUnstructuredGrid()  # Create a VTK unstructured grid

        print("vor_nodes", vor_nodes.shape)
        print("vor.vertices", vor.vertices.shape)
        print("vor.point_region", vor.point_region.shape)

        points_array = vtk.vtkPoints()  # Add the Voronoi vertices as points

        for v in vor.vertices:
            points_array.InsertNextPoint(v)
            # if crop_domain.dist(v) <= 0:
            #     # print(v)
            #     points_array.InsertNextPoint(v)
            # else:
            #     points_array.InsertNextPoint([0, 0, 0])  # simpler for vizualizing in ParaView

        # # print("\n additional ")
        # for n in add_nodes:  # additional nodes ...
        #     # print(n)
        #     points_array.InsertNextPoint(n)

        grid.SetPoints(points_array)

        vol = []
        for idx, reg_num in enumerate(vor.point_region):
            indices = vor.regions[reg_num]
            i = reg_num - 1
            if idx < noc:
                if not -1 in indices:  # not an open region
                    # inside = True
                    # for j in indices:
                    #     if crop_domain.dist(vor.vertices[j]) >= 0:  # negative values mean inside the domain (sdf)
                    #         inside = False
                    #         break
                    # if inside:
                    vol.append(ConvexHull(vor.vertices[indices]).volume)
                    id_array = vtk.vtkIdList()
                    for vertex_index in indices:
                        id_array.InsertNextId(vertex_index)
                    grid.InsertNextCell(vtk.VTK_CONVEX_POINT_SET, id_array)
                else:
                    vol.append(np.nan)
                    print("unbounded cell", i, ":", indices)

        cell_id = vtk.vtkDoubleArray()
        cell_id.SetName("cell_id")
        cell_id.SetNumberOfValues(noc)
        for j in range(0, noc):
            cell_id.SetValue(j, j)
        celldata = grid.GetCellData()
        celldata.AddArray(cell_id)

        vols = vtk.vtkDoubleArray()
        vols.SetName("volumes")
        vols.SetNumberOfValues(noc)
        for j in range(0, noc):
            vols.SetValue(j, vol[j])
        celldata = grid.GetCellData()
        celldata.AddArray(vols)

        return grid

    def make_periodic_(self, nodes, width_x, width_y):
        """maps the point periodically into [-width_x / 2., width_y / 2.] for x and y"""
        nodes_ = nodes.copy()
        nodes_[:, 0] = (nodes_[:, 0] + width_x / 2) % width_x - width_x / 2
        nodes_[:, 1] = (nodes_[:, 1] + width_y / 2) % width_y - width_y / 2
        return nodes_

    def get_cell_bounds(self, i: int, j: int, k: int):
        """returns the bounding box (minx ,miny ,minz, maxx, maxy, maxz) of soil cell with index i, j, k of a rectangular grid"""
        ms = self.ms  # rename
        width = ms.maxBound.minus(ms.minBound)
        cell_width = np.array([width.x / ms.resolution.x, width.y / ms.resolution.y, width.z / ms.resolution.z])
        min_ = np.array([ms.minBound.x + i * cell_width[0], ms.minBound.y + j * cell_width[1], ms.minBound.z + k * cell_width[2]])
        max_ = min_ + cell_width
        return min_, max_

    def nodes_within_box(self, nodes, min_, max_):
        """returns the nodes inside the bounding box given by min_, and max_ and the corresponding node indices ni"""
        mask = (
            (nodes[:, 0] >= min_[0])
            & (nodes[:, 0] < max_[0])
            & (nodes[:, 1] >= min_[1])
            & (nodes[:, 1] < max_[1])
            & (nodes[:, 2] >= min_[2])
            & (nodes[:, 2] < max_[2])
        )
        ni = np.where(mask)[0]
        return nodes[mask], ni

    def flip_(self, nodes, center, axis):
        """flips the nodes around the center according axis"""
        n_ = nodes - center
        n_[:, axis] = -n_[:, axis]
        n_ = n_ + center
        return n_

    def shift_(self, nodes, i: int, j: int, k: int):
        """adds shifted nodes to 8 sides horizontally plus two mirrored vertically with index i, j, k"""
        min_, max_ = self.get_cell_bounds(i, j, k)
        width_ = max_ - min_
        nodes_, ni = self.nodes_within_box(nodes, min_, max_)
        dx = width_[0]
        dy = width_[1]
        dz = width_[2]
        if nodes_.shape[0] > 0:
            flipped_z = self.flip_(nodes_, 0.5 * (min_ + max_), 2)  # flip them ...
            nodes_ = np.vstack((nodes_, nodes_ + np.array([dx, 0.0, 0.0])))
            nodes_ = np.vstack((nodes_, nodes_ + np.array([-dx, 0.0, 0.0])))
            nodes_ = np.vstack((nodes_, nodes_ + np.array([0.0, dy, 0.0])))
            nodes_ = np.vstack((nodes_, nodes_ + np.array([0.0, -dy, 0.0])))
            nodes_ = np.vstack((nodes_, nodes_ + np.array([dx, dy, 0.0])))
            nodes_ = np.vstack((nodes_, nodes_ + np.array([-dx, dy, 0.0])))
            nodes_ = np.vstack((nodes_, nodes_ + np.array([dx, -dy, 0.0])))
            nodes_ = np.vstack((nodes_, nodes_ + np.array([-dx, -dy, 0.0])))
            nodes_ = np.vstack((nodes_, flipped_z + np.array([0.0, 0.0, dz])))
            nodes_ = np.vstack((nodes_, flipped_z + np.array([0.0, 0.0, -dz])))
            return nodes_, ni
        else:  # no nodes within
            return np.ones((0,)), np.ones((0,))

    def mirror_(self, nodes, i: int, j: int, k: int):
        """adds mirrored nodes to the 6 sides of the cube with index i, j, k"""
        min_, max_ = self.get_cell_bounds(i, j, k)
        width_ = max_ - min_
        width_ = np.expand_dims(width_, axis=0)
        nodes_, ni = self.nodes_within_box(nodes, min_, max_)
        if nodes_.shape[0] > 0:
            fipped_n = [self.flip_(nodes_, 0.5 * (min_ + max_), i_) for i_ in range(0, 3)]  # flip them ...
            zeros = np.zeros((nodes_.shape[0], 1))
            translate_ = np.ones((nodes_.shape[0], 1)) @ width_
            trans0 = np.hstack((translate_[:, 0, np.newaxis], zeros, zeros))
            trans1 = np.hstack((zeros, translate_[:, 1, np.newaxis], zeros))
            trans2 = np.hstack((zeros, zeros, translate_[:, 2, np.newaxis]))
            nodes_ = np.vstack((nodes_, fipped_n[0] + trans0))  # add them
            nodes_ = np.vstack((nodes_, fipped_n[0] - trans0))
            nodes_ = np.vstack((nodes_, fipped_n[1] + trans1))
            nodes_ = np.vstack((nodes_, fipped_n[1] - trans1))
            nodes_ = np.vstack((nodes_, fipped_n[2] + trans2))
            nodes_ = np.vstack((nodes_, fipped_n[2] - trans2))
            return nodes_, ni
        else:  # no nodes within
            return np.ones((0,)), np.ones((0,))

    def get_default_volumes_(self):
        """returns the soil volumes of a rectangular grid (default)"""
        ms = self.ms  # rename
        width = ms.maxBound.minus(ms.minBound)
        vol = width.x * width.y * width.z / ms.resolution.x / ms.resolution.y / ms.resolution.z
        return np.ones((int(ms.resolution.x * ms.resolution.y * ms.resolution.z),)) * vol


if __name__ == "__main__":

    # Example A: create lookup table:

    # sand = [0.045, 0.43, 0.15, 3, 1000]
    # loam = [0.08, 0.43, 0.04, 1.6, 50]
    # clay = [0.1, 0.4, 0.01, 1.1, 10]
    hydrus_loam = [0.078, 0.43, 0.036, 1.56, 24.96]
    filename = "hydrus_loam"
    sp = vg.Parameters(hydrus_loam)
    vg.create_mfp_lookup(sp)
    peri = PerirhizalPython()
    #peri.create_lookup_mpi(filename, sp)  # takes some hours; mpiexec -n 4 python Perirhizal.py
    start = time.time()
    peri.create_lookup(filename, sp)  # should be better
    end = time.time()
    print("Creating the simplified lookup table for hydrus loam took ", end - start, " seconds.")
    start = time.time()
    peri.create_lookup_global(peri.water_filename)  # this is global
    end = time.time()
    print("Creating the general lookup table for all van Genuchten sets took ", end - start, " seconds.")
    
    # generate some tests
    

    # filename = "hydrus_loam"
    # sp = vg.Parameters(hydrus_loam)
    # vg.create_mfp_lookup(sp)
    # peri = PerirhizalPython()
    # peri.create_lookup_mpi(filename, sp)  # takes some hours; mpiexec -n 4 python Perirhizal.py

    # Example B: usage lookup table:
    # # peri.open_lookup(filename)
    # peri = PerirhizalPython()
    # loam = [0.08, 0.43, 0.04, 1.6, 50]
    # peri.set_soil(vg.Parameters(loam))
    # a = 0.1  # cm
    # kr = 1.73e-4  # [1/day]
    # rx = -15000  # cm
    # sx = -200  # cm
    # rho = 1 / a
    # inner_kr = a * kr
    # rsx = peri.soil_root_interface_potentials([rx], [sx], [inner_kr], [rho])
    # print("root soil interface", rsx, "cm")
    # print("results into a flux of", kr * 2 * a * np.pi * (rsx - rx), "cm3/day")

    # Example C: voronoi outer radii
    import matplotlib.pyplot as plt

    import plantbox.visualisation.vtk_plot as vp

    min_b = pb.Vector3d([-38.0, -8.0, -200.0])
    max_b = pb.Vector3d([38.0, 8.0, -0.301])  # <- exlude seed
    cell_number = pb.Vector3d([1, 1, 1])

    plant = pb.MappedPlant()
    path = "../../modelparameter/structural/plant/"
    plant.readParameters(path + "fspm2023" + ".xml")
    plant.initialize()
    plant.setRectangularGrid(min_b, max_b, cell_number, False, False)
    plant.simulate(40)
    # print(plant.nodes[0]) # seed is located at (0,0,-0.3)
    # vp.plot_plant(plant, "subType")

    peri = PerirhizalPython(plant)
    outer_radii = peri.get_outer_radii("voronoi_periodic")  # length, surface, volume, voronoi_periodic, voronoi_bounded
    print("outer_radii.shape:", outer_radii.shape, "np.nanmin(outer_radii):", np.nanmin(outer_radii), "np.nanmax(outer_radii):", np.nanmax(outer_radii))

    lengths = np.array(plant.segLength())
    print("lengths.shape:", len(lengths), "np.nanmin(lengths):", np.nanmin(lengths), "np.nanmax(lengths):", np.nanmax(lengths))
    inner_radii = np.array(plant.getEffectiveRadii())
    inner_radii[outer_radii == 0] = 0  # if outer radius is zero, there is no perirhizal zone, and the inner radius must be set to zero as well
    vol = np.prod(plant.getDomainWidth())
    print("Domain volume", vol, "cm3")
    print("total perirhizal volume", np.pi * np.sum(outer_radii * outer_radii * lengths) - np.pi * np.sum(inner_radii * inner_radii * lengths), "cm3")
    # if the seed is inside the domain voronoi is not exact due to the seed node voronoi volume (which apical to any root semgent)

    outer_radii = np.clip(outer_radii, 0, 2)  # to avoid outliers for better vizualization
    plt.hist(outer_radii, weights=lengths, bins=40, rwidth=0.9)
    plt.show()
