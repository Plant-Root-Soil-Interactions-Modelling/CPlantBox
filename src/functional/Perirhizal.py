import numpy as np
from numpy import linalg as LA
import scipy.linalg as la
from scipy.interpolate import RegularGridInterpolator
from scipy.optimize import fsolve, root_scalar
from scipy.spatial import ConvexHull, Voronoi

import plantbox as pb
import plantbox.functional.van_genuchten as vg
from plantbox import MappedSegments, Perirhizal


class PerirhizalPython(Perirhizal):
    """
    Helper class for modelling the perirhizal zone

    Wraps MappedSegments (or specialisations MappedPlant) and adds functions to retrieve information

    * calculates root soil interface potential using steady rate approximation (Schröder et al. 2008)
    * perirhizal_conductance_per_layer calculates the perirhizal conductance in a soil layer (or cell) (Vanderborght et al. 2023, Eqn [6])
    * support of 4D lookup tables (file type is a zipped archive)
    * analysis across soil grid, e.g. get_density, average, aggregate
    * calculates outer perirhizal radii (based on denisties or Voronoi)

    run script to create a 4D lookup table (see __main__)
    """

    def __init__(self, ms=None):
        """ms      reference to MappedSegments"""
        if ms:
            super().__init__(ms)
        else:
            super().__init__()

        self.lookup_table = None  # optional 4d look up table to find soil root interface potentials
        self.simp_lookup_table = None  # optional 2d look up table to find soil root interface potentials
        self.global_lookup_table = None  # optional 3d look up table to find soil root interface potentials
        self.sp = None  # corresponding van gencuchten soil parameter
        
        self.alpha_0 = 0.3 # a constant that is used for numerically solving the perirhizal waterflow
        self.h_wilt = -16000
        self.water_filename = "lookup_perirhizal_waterflow_global" #name for the global lookup table

    def set_soil(self, sp):
        """sets VG parameters, and no look up table (slow)"""
        vg.create_mfp_lookup(sp)
        self.sp = sp
        self.lookup_table = None
        self.simp_lookup_table = None

    def open_lookup(self, filename):
        """  opens a look-up table from a file to quickly find soil root interface potentials  """
        npzfile = np.load(filename + ".npz")
        interface = npzfile["interface"]
        rx_, sx_, akr_, rho_ = npzfile["rx_"], npzfile["sx_"], npzfile["akr_"], npzfile["rho_"]
        soil = npzfile["soil"]
        self.lookup_table = RegularGridInterpolator((rx_, sx_, akr_, rho_), interface)  # method = "nearest" fill_value = None , bounds_error=False
        self.sp = vg.Parameters(soil)
        
    def open_simp_lookup(self, filename):
        """opens a smaller look-up table from a file to quickly find soil root interface potentials"""
        npzfile = np.load(filename + "_simp.npz")
        interface = npzfile["interface"]
        inner_kr_b, base_mfp = npzfile["inner_kr_b_"], npzfile["base_mfp_"]
        print(inner_kr_b[0],inner_kr_b[-1], base_mfp[0],base_mfp[-1])
        soil = npzfile["soil"]
        self.simp_lookup_table = RegularGridInterpolator((inner_kr_b, base_mfp), interface)  # method = "nearest" fill_value = None , bounds_error=False
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
        
        rsx = np.array([PerirhizalPython.soil_root_interface_(rx[i], sx[i], inner_kr[i], rho[i], self.sp) for i in range(0, len(rx))])
        print("Norm of the computed interface potentials:", LA.norm(rsx))
        
        rsx1 = 0*rsx.copy()
        rsx2 = 0*rsx.copy()
        rsx3 = 0*rsx.copy()
        
        
        if self.lookup_table:
            rsx1 = self.soil_root_interface_potentials_table(rx, sx, inner_kr, rho)
            print("Norm of the difference basic lookup table:", LA.norm(rsx-rsx1))
        if self.simp_lookup_table:
            rsx2 = self.soil_root_interface_potentials_simp_table(rx, sx, inner_kr, rho)
            print("Norm of the difference simp lookup table:", LA.norm(rsx-rsx2))
        if self.global_lookup_table:
            rsx3 = self.soil_root_interface_potentials_global_table(rx, sx, inner_kr, rho)
            print("Norm of the difference global lookup table:", LA.norm(rsx-rsx3))
        
        return rsx, rsx1, rsx2, rsx3

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
        k_soilfun = lambda hsoil, hint: (vg.fast_mfp[sp](hsoil) - vg.fast_mfp[sp](hint)) / (hsoil - hint + 0.001)  # Vanderborgth et al. 2023, Eqn [7]
        rho2 = np.square(rho)  # rho squared
        b = 2 * (rho2 - 1) / (1 - 0.53 * 0.53 * rho2 + 2 * rho2 * (np.log(rho) + np.log(0.53)))  # Vanderborgth et al. 2023, Eqn [8]
        fun = lambda x: (inner_kr * rx + b * sx * k_soilfun(sx, x)) / (b * k_soilfun(sx, x) + inner_kr) - x
        #print(fun(rx),fun(sx-0.1),fun(sx+0.1), rho)
        if rx == sx:  # degenerate bracket: no gradient, interface equals bulk potential
            return sx  
        rsx = root_scalar(fun, method="brentq", bracket=[min(rx, sx), max(rx, sx)])
        return rsx.root
     
    @staticmethod
    def soil_root_interface_simp_(self, inner_kr_b, base_mfp, sp):
        """
        finds matric potential at the soil root interface for as single segment

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
        print(fun(x_int[0]),fun(x_int[1]), inner_kr_b, base_mfp)
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
        
        (sp_base is the van genuchten parameter set of Ks = 1 and alpha = alpha_0, m is an input)
        """
        
        
        
        fun = lambda x: self.lookup_global_mfp((vg_m, x)) + inner_kr_b * x +base_mfp
        x_int = [-15990.0,-0.1] # interval for search, get closer to 0 as x = h_sr * alpha / alpha_0, alpha / alpha_0 <= 1
        if fun(x_int[0])> 0:
            return x_int[0]
        if fun(x_int[1])< 0:
            return x_int[1]
        rsx = root_scalar(fun, method="brentq", bracket=[x_int[0], x_int[1]])
        #except:
        #    print("f_values", fun(x_int[0]), fun(x_int[1]), vg_m)
        #    rsx = 0
        return rsx.root
        

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

    def create_lookup_mpi(self, filename, sp):
        """
        Precomputes all soil root interface potentials for a specific soil type
        and saves results into a 4D lookup table

        filename       three files are written (filename, filename_, and filename_soil)
        sp             van genuchten soil parameters, , call
                       vg.create_mfp_lookup(sp) before
        """
        import os

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
        rho_ = np.logspace(np.log10(2.), np.log10(200.), rhon) # TODO: talk about the lower bound here: rho = 1 leads to problems
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
    
    def create_lookup_simp(self, filename, sp):
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
        
        print(inner_kr_b_int)
        min_int, max_int = vg.fast_mfp[sp](-sx_int_abs[1]), vg.fast_mfp[sp](-sx_int_abs[0])
        print(vg.fast_mfp[sp](-sx_int_abs[1]))
        print(vg.fast_mfp[sp](-sx_int_abs[0]))
        base_mfp_n_1 = 500
        base_mfp_n_2 = 500
        base_mfp_n = base_mfp_n_1 + base_mfp_n_2
        base_mfp_int = [inner_kr_b_int[0]*(0*self.h_wilt-(-rx_int_abs[0])) - vg.fast_mfp[sp](-sx_int_abs[0]), inner_kr_b_int[1]*(0*self.h_wilt-(-rx_int_abs[1])) - vg.fast_mfp[sp](-sx_int_abs[1])]
        #base_mfp_int are negative
        #base_mfp_int[1] = min(base_mfp_int[1],-1.0e-9) #ensure positivity, this can be close to 0.0, so do not take a huge tolerance here
        base_mfp_interval1= -np.logspace(np.log10(-base_mfp_int[0]), np.log10(1.0e-9), base_mfp_n_1)
        base_mfp_interval2= np.logspace(np.log10(1.0e-9), np.log10(base_mfp_int[1]), base_mfp_n_2)
        print(base_mfp_interval1,base_mfp_interval2)
        base_mfp_ = np.concatenate((base_mfp_interval1,base_mfp_interval2))
        #base_mfp_= np.linspace(base_mfp_int[0], base_mfp_int[1], base_mfp_n)
        print(base_mfp_int)
        interface = np.zeros((inner_kr_bn,base_mfp_n))
        
        tol = 0.001 #tolerance for extreme values of the matrix flux potential at the soil -> output h_rs = h_soil or h_x
        
        max_int = vg.fast_mfp[sp](-1)
        min_int = vg.fast_mfp[sp](-16000+1)
        
        for i, inner_kr_b in enumerate(inner_kr_b_):
            print(i)
            for j, base_mfp in enumerate(base_mfp_):
                #print(base_mfp, tol - (-1-self.h_wilt)*inner_kr_b - max_int, - tol - (-16000-self.h_wilt)*inner_kr_b - min_int, min_int)
                #base_mfp = max(base_mfp, tol - (-1-0*self.h_wilt)*inner_kr_b - max_int)
                #base_mfp = min(base_mfp, - tol - (-16000-0*self.h_wilt)*inner_kr_b - min_int)
                #print(base_mfp)
                interface[i, j] = PerirhizalPython.soil_root_interface_simp_(self, inner_kr_b, base_mfp, sp)
        np.savez(filename+"_simp", interface=interface, inner_kr_b_=inner_kr_b_, base_mfp_=base_mfp_, soil=list(sp))
        print(inner_kr_b_)
        print(base_mfp_)
        self.simp_lookup_table = RegularGridInterpolator((inner_kr_b_, base_mfp_), interface)
        self.sp = sp

    def create_lookup_global(self, filename, sp):
        """
        Precomputes all soil root interface potentials for a specific soil type
        and saves results into a 3D lookup table

        filename       three files are written (filename, filename_, and filename_soil)
        sp             van genuchten soil parameters, , call
                       vg.create_mfp_lookup(sp) before
        """

        
        # intervals for all inputs
        rx_int_abs = [1.0, 16000.0] #absolute values for logarithmic scaling
        sx_int_abs = [1.0, 16000.0]
        akr_int = [1.0e-7,1.0e-4]
        rho_int = [2.0,200.0]
        vg_m_int = [0.1,0.8]
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
        
        #construct a smaller lookup table for the simplified 
        for i, vg_m in enumerate(vg_m_):
            print(i)
            sp_dummy = vg.Parameters([0.1,0.4,self.alpha_0,1./(1.-vg_m),1.0])
            vg.create_mfp_lookup(sp_dummy)
            for j, sx in enumerate(sx_):
                #interface_vg[i, j] = vg.fast_mfp[sp_dummy](sx)
                interface_vg[i, j] = vg.matric_flux_potential(sx, sp_dummy)
            print(sx, vg.matric_flux_potential(sx+15000, sp_dummy), vg.matric_flux_potential(-1, sp_dummy))
        #np.savez(filename, interface=interface_vg, vg_m_=vg_m_, sx_=sx_, soil=list(sp))
        self.lookup_global_mfp = RegularGridInterpolator((vg_m_, sx_), interface_vg)        
        print(vg_m_,sx_)
        inner_kr_bn = 200
        inner_kr_b_abs_int = [akr_int[0]/(b_int[1]*vg_Ks_int[1]),akr_int[1]/(b_int[0]*vg_Ks_int[0])]
        inner_kr_b_= np.logspace(np.log10(inner_kr_b_abs_int[0]), np.log10(inner_kr_b_abs_int[1]), inner_kr_bn)
        print(inner_kr_b_abs_int)
        
        #base_mfp_n = 200
        base_mfp_n_1 = 100
        base_mfp_n_2 = 100
        base_mfp_n = base_mfp_n_1 + base_mfp_n_2
        print(-sx_int_abs[1] * vg_alpha_int[1]/self.alpha_0)
        print(self.lookup_global_mfp((vg_m_int[0],-sx_int_abs[1] * vg_alpha_int[1]/self.alpha_0*0.1)))
        print(vg.fast_mfp[sp_dummy](-10000.0))
        print("obere grenze",-sx_int_abs[0] * vg_alpha_int[0]/self.alpha_0)
        print(self.lookup_global_mfp((vg_m_int[1],min(-sx_int_abs[0] * vg_alpha_int[0]/self.alpha_0,-1.0))))
        base_mfp_int = [inner_kr_b_abs_int[0]*(-(-rx_int_abs[0]) * vg_alpha_int[0]/self.alpha_0) - self.lookup_global_mfp((vg_m_int[1],min(-sx_int_abs[0] * vg_alpha_int[0]/self.alpha_0,-1))), inner_kr_b_abs_int[1]*(-(-rx_int_abs[1]) *vg_alpha_int[1]/self.alpha_0) -self.lookup_global_mfp((vg_m_int[0],-sx_int_abs[1] * vg_alpha_int[1]/self.alpha_0))] 
        base_mfp_interval1= -np.logspace(np.log10(-base_mfp_int[0]), np.log10(1.0e-9), base_mfp_n_1)
        base_mfp_interval2= np.logspace(np.log10(1.0e-9), np.log10(base_mfp_int[1]), base_mfp_n_2)
        print(base_mfp_interval1,base_mfp_interval2)
        base_mfp_ = np.concatenate((base_mfp_interval1,base_mfp_interval2))
        #print(base_mfp_abs_int)
        #base_mfp_= - np.logspace(np.log10(-base_mfp_abs_int[0]), np.log10(-base_mfp_abs_int[1]), base_mfp_n)
        
        interface = np.zeros((vg_m_n,inner_kr_bn,base_mfp_n))
        
        tol = 0.001 #tolerance for extreme values of the matrix flux potential at the soil -> output h_rs = h_soil or h_x
        
        
        for i, vg_m in enumerate(vg_m_):
            print(i)
            max_int = self.lookup_global_mfp((vg_m,-1))
            min_int = self.lookup_global_mfp((vg_m,-16000))
            for j, inner_kr_b in enumerate(inner_kr_b_):
                for k, base_mfp in enumerate(base_mfp_):
                    #base_mfp = max(base_mfp, tol - (-1)*inner_kr_b - max_int)
                    #base_mfp = min(base_mfp, -tol - (-16000)*inner_kr_b - min_int)
                    interface[i, j, k] = PerirhizalPython.soil_root_interface_global_(self, inner_kr_b, base_mfp, vg_m)
                    #if base_mfp + (-1)*inner_kr_b + max_int < tol:
                    #    interface[i, j, k] = -1
                    #else:
                    #    if base_mfp + (-16000)*inner_kr_b + min_int > tol:
                    #        interface[i, j, k] = -16000
                    #    else:
                    #        interface[i, j, k] = PerirhizalPython.soil_root_interface_global_(self, inner_kr_b, base_mfp, vg_m)
        np.savez(filename, interface=interface, interface_vg=interface_vg, vg_m_=vg_m_, inner_kr_b_=inner_kr_b_, base_mfp_=base_mfp_, sx_=sx_, soil=list(sp))
        print(vg_m_, inner_kr_b_, base_mfp_)
        self.global_lookup_table = RegularGridInterpolator((vg_m_, inner_kr_b_, base_mfp_), interface)
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
        try:
            sx = np.array(sx)
            mask = inner_kr_ == 0
            inner_kr_[mask] = 1.0e-7
            rsx = self.lookup_table((rx, sx, inner_kr_, rho_))
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

    def soil_root_interface_potentials_simp_table(self, rx, sx, inner_kr_, rho_):
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
            print(b)
            inner_kr_b = np.divide(inner_kr_,b)
            print(inner_kr_b)
            base_mfp = - (inner_kr_b * rx + vg.fast_mfp[self.sp](sx))
            print(vg.fast_mfp[self.sp](sx))
            print(base_mfp)
            
            #tol = 0.01
            #mask2 = base_mfp - tol > 0
            #base_mfp[mask2] = -tol
            
            rsx = self.simp_lookup_table((inner_kr_b, base_mfp))
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
        
    def soil_root_interface_potentials_global_table(self, rx, sx, inner_kr_, rho_):
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
            #print("inner_kr1",inner_kr_[0]*self.sp.alpha/self.sp.Ksat/b[0]/self.alpha_0)
            inner_kr_b = np.divide(inner_kr_,self.sp.Ksat*b)
            base_mfp = - ( inner_kr_b * rx * self.sp.alpha / self.alpha_0 + self.lookup_global_mfp((self.sp.m,sx*self.sp.alpha/self.alpha_0)))
            print((inner_kr_b, base_mfp, self.sp.m)) #TODO multiply by alpha or something
            rsx = np.array([self.global_lookup_table((self.sp.m, inner_kr_b[i], base_mfp[i])) for i in range(len(inner_kr_b))])*self.alpha_0/self.sp.alpha
            print(mask, rsx)
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
            return self.get_outer_radii_voronoi()  # default is bounded voronoi
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
        ms = self.ms  # renamethat
        cell2seg = ms.cell2seg
        radii = ms.getEffectiveRadii()
        lengths = ms.segLength()
        width = ms.maxBound.minus(ms.minBound)
        outer_r = np.zeros((len(radii),))
        if volumes is None:
            n = ms.resolution.x * ms.resolution.y * ms.resolution.z
            volumes = np.ones(int(n)) * (width.x * width.y * width.z) / n
            # print("PerirhizalPython.get_outer_radii_(): each soil volume has", (width.x * width.y * width.z) / n, "cm3")

        for cell_id, seg_ids in cell2seg.items():
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
        min_b = ms.minBound
        max_b = ms.maxBound

        ns = ms.getNumberOfMappedSegments()
        nodes_ = ms.nodes
        nodes = np.array([[n.x, n.y, n.z] for n in nodes_])  # to numpy array
        nodes = nodes[: ns + 1]  # TODO there are more nodes than segments,
        # e.g. some nodes are not connected to the root system.For now we just drop them, but we should check why they are there in the first place.
        print("nodes", len(nodes_))
        print("nodes truncated", nodes.shape)
        print("segments", ns)

        x = int(ms.resolution.x)
        y = int(ms.resolution.y)
        z = int(ms.resolution.z)

        width = np.array([max_b.x, max_b.y, max_b.z]) - np.array([min_b.x, min_b.y, min_b.z])
        # print("making periodic", width)
        # print("nodes", nodes.shape)
        # print("resolution", x, y, z)
        nodes = self.make_periodic_(nodes, width[0], width[1])

        vol = np.empty((nodes.shape[0]))
        vol[:] = np.nan

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

    def get_voronoi_mesh(self):
        """Creates a VTK grid containing vor_nodes; assumes periodic boundaries in x and y direction, and open in z direction"""
        nodes_ = self.ms.nodes  # make nodes periodic and pass all
        nodes = np.array([[n.x, n.y, n.z] for n in nodes_])  # to numpy array
        min_b = self.ms.minBound
        max_b = self.ms.maxBound
        width = np.array([max_b.x, max_b.y, max_b.z]) - np.array([min_b.x, min_b.y, min_b.z])
        nodes = self.make_periodic_(nodes, width[0], width[1])
        return self.get_voronoi_mesh_(nodes, len(nodes))

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
        """adds shifted nodes to the 6 sides of the cube with index i, j, k"""
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

    sand = [0.045, 0.43, 0.15, 3, 1000]
    loam = [0.08, 0.43, 0.04, 1.6, 50]
    clay = [0.1, 0.4, 0.01, 1.1, 10]

    hydrus_loam = [0.078, 0.43, 0.036, 1.56, 24.96]
    hydrus_clay = [0.068, 0.38, 0.008, 1.09, 4.8]
    hydrus_sand = [0.045, 0.43, 0.145, 2.68, 712.8]
    hydrus_sandyloam = [0.065, 0.41, 0.075, 1.89, 106.1]

    filename = "hydrus_loam"
    sp = vg.Parameters(hydrus_loam)
    vg.create_mfp_lookup(sp)
    peri = PerirhizalPython()
    #peri.create_lookup_mpi(filename, sp)  # takes some hours; mpiexec -n 4 python Perirhizal.py
    peri.create_lookup_simp(filename, sp)  # should be better
    peri.create_lookup_global(peri.water_filename)  # this is global
    
    # generate some tests
    

    # # peri.open_lookup(filename)
    # peri.set_soil(vg.Parameters(loam))
    # a = 0.1  # cm
    # kr = 1.73e-4  # [1/day]
    # rx = -15000  # cm
    # sx = 0.0  # cm
    # rho = 1 / a
    # inner_kr = a * kr
    # rsx = peri.soil_root_interface_potentials([rx], [sx], [inner_kr], [rho])
    # print("root soil interface", rsx, "cm")
    # print("results into a flux of", kr * 2 * a * np.pi * (rsx - rx), "cm3/day")
