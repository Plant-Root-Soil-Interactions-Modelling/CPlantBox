import sys; sys.path.append(".."); sys.path.append("../..");

from plantbox import Perirhizal
import functional.van_genuchten as vg

import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.optimize import fsolve


class PerirhizalPython(Perirhizal):
    """
    Helper class for modelling perirhizal zone 
    
    * calculates outer perirhizal radii    
    * calculates root soil interface potential using steady rate approximation (Schröder et al. 2008)
    * support of 4d lookup tables (file type is a zipped archive)
    
    """

    def __init__(self, ms = None):
        """ 
        ms      reference to MappedSegments
        """
        if ms:
            super().__init__(ms)
        else:
            super().__init__()

        self.lookup_table = None  # optional 4d look up table to find soil root interface potentials
        self.sp = None  # corresponding van gencuchten soil parameter

    def seg_outer_radii(self, type:int, vols = None):
        """
        calculates the outer perirhizal radii assuming homogeneously distributed roots within soil volumes,
        calls Perirhizal.cpp
        
        type             assigns the radii proportional to segment volume (0) segment surface (1), or segment length (2)
        vols             (optional) cell volumes in case of not equidistant grids         
        """
        if vols:
            super().segOuterRadii(type, vols)
        else:
            super().segOuterRadii(type)

    def soil_root_interface_potentials(self, rx, sx, inner_kr, rho, sp):
        """
        finds matric potentials at the soil root interface for as all segments
        
        rx             xylem matric potential [cm]
        sx             bulk soil matric potential [cm]
        inner_kr       root radius times hydraulic conductivity [cm/day] 
        rho            geometry factor [1] (outer_radius / inner_radius)
        sp             soil van Genuchten parameters (type vg.Parameters), call 
                       vg.create_mfp_lookup(sp) before 
        """
        assert len(rx) == len(sx) == len(inner_kr) == len(rho), "rx, sx, inner_kr, and rho must have the same length"
        if self.lookup_table:
            rsx = self.soil_root_interface_potentials_table(rx, sx, inner_kr_, rho)
        else:
            rsx = np.array([PerirhizalPython.soil_root_interface_(rx[i], sx[i], inner_kr[i], rho[i], sp) for i in range(0, len(rx))])
        return rsx

    @staticmethod
    def soil_root_interface_(rx, sx, inner_kr, rho, sp):
        """
        finds matric potential at the soil root interface for as single segment
        
        rx             xylem matric potential [cm]
        sx             bulk soil matric potential [cm]
        inner_kr       root radius times hydraulic conductivity [cm/day] 
        rho            geometry factor [1] (outer_radius / inner_radius)
        sp             soil van Genuchten parameters (type vg.Parameters), call 
                       vg.create_mfp_lookup(sp) before 
        """
        k_soilfun = lambda hsoil, hint: (vg.fast_mfp[sp](hsoil) - vg.fast_mfp[sp](hint)) / (hsoil - hint)
        # rho = outer_r / inner_r  # Eqn [5]
        rho2 = rho * rho  # rho squaredswitched
        # b = 2 * (rho2 - 1) / (1 + 2 * rho2 * (np.log(rho) - 0.5))  # Eqn [4]
        b = 2 * (rho2 - 1) / (1 - 0.53 * 0.53 * rho2 + 2 * rho2 * (np.log(rho) + np.log(0.53)))  # Eqn [7]
        fun = lambda x: (inner_kr * rx + b * sx * k_soilfun(sx, x)) / (b * k_soilfun(sx, x) + inner_kr) - x
        rsx = fsolve(fun, (rx + sx) / 2)
        return rsx

    def create_lookup(self, filename, sp):
        """      
        Precomputes all soil root interface potentials for a specific soil type 
        and saves results into a 4D lookup table
        
        filename       three files are written (filename, filename_, and filename_soil)
        sp             van genuchten soil parameters, , call 
                       vg.create_mfp_lookup(sp) before 
        """
        rxn = 150
        rx_ = -np.logspace(np.log10(1.), np.log10(16000), rxn)
        rx_ = rx_ + np.ones((rxn,))
        rx_ = rx_[::-1]
        sxn = 150
        sx_ = -np.logspace(np.log10(1.), np.log10(16000), sxn)
        sx_ = sx_ + np.ones((sxn,))
        sx_ = sx_[::-1]
        akrn = 100
        akr_ = np.logspace(np.log10(1.e-7), np.log10(1.e-4), akrn)
        rhon = 30
        rho_ = np.logspace(np.log10(1.), np.log10(200.), rhon)
        interface = np.zeros((rxn, sxn, akrn, rhon))
        for i, rx in enumerate(rx_):
            print(i)
            for j, sx in enumerate(sx_):
                for k, akr in enumerate(akr_):
                    for l, rho in enumerate(rho_):
                        interface[i, j, k, l] = PerirhizalPython.soil_root_interface_(rx, sx, akr, rho, sp)
        np.savez(filename, interface = interface, rx_ = rx_, sx_ = sx_, akr_ = akr_, rho_ = rho_, soil = list(sp))
        self.lookup_table = RegularGridInterpolator((rx_, sx_, akr_, rho_), interface)
        self.sp = sp

    def open_lookup(self, filename):
        """ 
        Opens a look-up table from a file to quickly find soil root interface potentials 
        """
        npzfile = np.load(filename + ".npz")
        interface = npzfile["interface"]
        rx_, sx_, akr_, rho_ = npzfile["rx_"], npzfile["sx_"], npzfile["akr_"], npzfile["rho_"]
        soil = npzfile["soil"]
        self.lookup_table = RegularGridInterpolator((rx_, sx_, akr_, rho_), interface)  # method = "nearest" fill_value = None , bounds_error=False
        self.sp = vg.Parameters(soil)

    def soil_root_interface_potentials_table(self, rx, sx, inner_kr_, rho_):
        """
        finds potential at the soil root interface using a lookup table
            
        rx             xylem matric potential [cm]
        sx             bulk soil matric potential [cm]
        inner_kr       root radius times hydraulic conductivity [cm/day] 
        rho            geometry factor [1]
        f              function to look up the potentials
        """
        try:
            rsx = self.lookup_table((rx, sx, inner_kr_ , rho_))
        except:
            print("Look up failed: ")
            print("rx", np.min(rx), np.max(rx))  # 0, -16000
            print("sx", np.min(sx), np.max(loamsx))  # 0, -16000
            print("inner_kr", np.min(inner_kr_), np.max(inner_kr_))  # 1.e-7 - 1.e-4
            print("rho", np.min(rho_), np.max(rho_))  # 1. - 200.
        return rsx


if __name__ == "__main__":

    sand = [0.045, 0.43, 0.15, 3, 1000]
    loam = [0.08, 0.43, 0.04, 1.6, 50]
    clay = [0.1, 0.4, 0.01, 1.1, 10]
    hydrus_loam = [0.078, 0.43, 0.036, 1.56, 24.96]
    hydrus_clay = [0.068, 0.38, 0.008, 1.09, 4.8]
    hydrus_sand = [0.045, 0.43, 0.145, 2.68, 712.8]
    hydrus_sandyloam = [0.065, 0.41, 0.075, 1.89, 106.1]

    # e.g. for loam
    filename = "table_loam"
    sp = vg.Parameters(loam)
    vg.create_mfp_lookup(sp)
    peri = PerirhizalPython()
    peri.create_lookup(filename, sp)  # takes some hours

    peri.open_lookup(filename)
