import sys; sys.path.append(".."); sys.path.append("../..");

import plantbox as pb
from plantbox import Perirhizal
from plantbox import MappedSegments
import functional.van_genuchten as vg

import numpy as np

from scipy.spatial import ConvexHull, Voronoi
from scipy.interpolate import RegularGridInterpolator
from scipy.optimize import fsolve
from scipy.optimize import root_scalar
import scipy.linalg as la


class PerirhizalPython(Perirhizal):
    """
    Helper class for modelling the perirhizal zone 
    
    Wraps MappedSegments (or specialisations MappedPlant) and adds functions to retrieve information 
    
    * calculates outer perirhizal radii    
    * calculates root soil interface potential using steady rate approximation (Schröder et al. 2008)
    * support of 4d lookup tables (file type is a zipped archive)
    
    run script to create a 4d lookup table (see __main__)
    """

    def __init__(self, ms = None):
        """  ms      reference to MappedSegments """
        if ms:
            super().__init__(ms)
        else:
            super().__init__()

        self.lookup_table = None  # optional 4d look up table to find soil root interface potentials
        self.sp = None  # corresponding van gencuchten soil parameter

    def set_soil(self, sp):
        """ sets VG parameters, and no look up table (slow) """
        vg.create_mfp_lookup(sp)
        self.sp = sp
        self.lookup_table = None

    def open_lookup(self, filename):
        """  opens a look-up table from a file to quickly find soil root interface potentials  """
        npzfile = np.load(filename + ".npz")
        interface = npzfile["interface"]
        rx_, sx_, akr_, rho_ = npzfile["rx_"], npzfile["sx_"], npzfile["akr_"], npzfile["rho_"]
        soil = npzfile["soil"]
        self.lookup_table = RegularGridInterpolator((rx_, sx_, akr_, rho_), interface)  # method = "nearest" fill_value = None , bounds_error=False
        self.sp = vg.Parameters(soil)

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
        if self.lookup_table:
            rsx = self.soil_root_interface_potentials_table(rx, sx, inner_kr, rho)
        else:
            rsx = np.array([PerirhizalPython.soil_root_interface_(rx[i], sx[i], inner_kr[i], rho[i], self.sp) for i in range(0, len(rx))])
            # rsx = rsx[:, 0]
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
        if inner_kr < 1.e-7:
            return sx
        k_soilfun = lambda hsoil, hint: (vg.fast_mfp[sp](hsoil) - vg.fast_mfp[sp](hint)) / (hsoil - hint)  # Vanderborgth et al. 2023, Eqn [7]
        rho2 = np.square(rho)  # rho squared
        b = 2 * (rho2 - 1) / (1 - 0.53 * 0.53 * rho2 + 2 * rho2 * (np.log(rho) + np.log(0.53)))  # Vanderborgth et al. 2023, Eqn [8]
        fun = lambda x: (inner_kr * rx + b * sx * k_soilfun(sx, x)) / (b * k_soilfun(sx, x) + inner_kr) - x
        rsx = root_scalar(fun, method = 'brentq', bracket = [min(rx, sx), max(rx, sx)])
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
        r_phriz = np.divide(1., np.sqrt(np.pi * rld))  # outer perirhizal radii per layer [cm]
        rho = np.divide(r_phriz, r_root)  # [1], see Vanderborght et al. (2023), Eqn [9]
        rho2 = np.square(rho)  # [1]
        b = np.divide(2 * (rho2 - np.ones(rho2.shape)) , np.ones(rho2.shape) - 0.53 * 0.53 * rho2 + 2 * rho2 * (np.log(rho) + np.log(0.53)))  # [1], see Eqn [8]
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
        from mpi4py import MPI
        import os
        comm = MPI.COMM_WORLD
        size = comm.Get_size()
        rank = comm.Get_rank()

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
                print('at index', index_ + 1 , "/", stop - start, "on thread", rank)

            rx = rx_[i]
            sx = sx_[j]
            akr = akr_[k]
            rho = rho_[l]

            interface_local[index_] = PerirhizalPython.soil_root_interface_(rx, sx, akr, rho, sp)

        # data2share needs to use floats for Allgatherv with MPI.DOUBLE to work.
        data2share = np.array(interface_local, dtype = np.float64)

        # other data needed by comm.Allgatherv
        all_sizes = np.full(size, count)
        all_sizes[:remainder] += 1

        offsets = np.zeros(len(all_sizes), dtype = np.int64)
        offsets[1:] = np.cumsum(all_sizes)[:-1]
        all_sizes = tuple(all_sizes)
        offsets = tuple(offsets)

        # share the vectors
        interface = np.zeros(work_size)
        comm.Allgatherv([data2share, MPI.DOUBLE], [interface, all_sizes, offsets, MPI.DOUBLE])

        if rank == 0:
            interface = interface.reshape(rxn, sxn, akrn, rhon)  # reset shape
            np.savez(filename, interface = interface, rx_ = rx_, sx_ = sx_, akr_ = akr_, rho_ = rho_, soil = list(sp))
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

            if np.max(rx) > 0: print("xylem matric potential positive", np.max(rx), "at", np.argmax(rx))
            if np.min(rx) < -16000: print("xylem matric potential under -16000 cm", np.min(rx), "at", np.argmin(rx))
            if np.max(sx) > 0: print("soil matric potential positive", np.max(sx), "at", np.argmax(sx))
            if np.min(sx) < -16000: print("soil matric potential under -16000 cm", np.min(sx), "at", np.argmin(sx))
            if np.min(inner_kr_) < 1.e-7: return sx  # for kr ~ 0, rsx == sx
            if np.max(inner_kr_) > 1.e-4: print("radius times radial conductivity above 1.e-4", np.max(inner_kr_), "at", np.argmax(inner_kr_))
            if np.min(rho_) < 1: print("geometry factor below 1", np.min(rho_), "at", np.argmin(rho_))
            if np.max(rho_) > 200: print("geometry factor above 200", np.max(rho_), "at", np.argmax(rho_))

            print("PerirhizalPython.soil_root_interface_potentials_table(): table look up failed, value exceeds table")
            print("\trx", np.min(rx), np.max(rx), "sx", np.min(sx), np.max(sx), "inner_kr", np.min(inner_kr_), np.max(inner_kr_), "rho", np.min(rho_), np.max(rho_))

        return rsx

    def get_cell_bounds(self, i:int, j:int, k:int):
        """ returns the bounding box (minx ,miny ,minz, maxx, maxy, maxz) of soil cell with index i, j, k of a rectangular grid """
        ms = self.ms  # rename
        width = ms.maxBound.minus(ms.minBound)
        cell_width = np.array([width.x / ms.resolution.x, width.y / ms.resolution.y, width.z / ms.resolution.z])
        min_ = np.array([ms.minBound.x + i * cell_width[0], ms.minBound.y + j * cell_width[1], ms.minBound.z + k * cell_width[2]])
        max_ = min_ + cell_width
        return min_, max_

    def nodes_within_box(self, nodes, min_, max_):
        """ returns the nodes inside the bounding box given by min_, and max_ """
        one = np.ones((nodes.shape[0], 1))
        n_, ni = [], []
        for i, n in enumerate(nodes):
            if n[0] >= min_[0] and n[1] >= min_[1] and n[2] >= min_[2] and n[0] < max_[0] and n[1] < max_[1] and n[2] < max_[2]:
                n_.append(n)
                ni.append(i)
        # print("nodes within", min_, max_, ":", len(n_))
        return np.array(n_), np.array(ni)

    def mirror(self, nodes, i:int, j:int, k:int):
        """ adds mirrored nodes to the 6 sides of the cube with index i, j, k"""
        min_, max_ = self.get_cell_bounds(i, j, k)
        width_ = max_ - min_
        width_ = np.expand_dims(width_, axis = 0)
        center_ = min_ + width_ / 2.
        n, ni = self.nodes_within_box(nodes, min_, max_)
        if n.shape[0] > 0:
            nodes_surr = n
            fipped_n = [self.flip_(n, center_, i_) for i_ in range(0, 3)]  # flip them ...
            zeros = np.zeros((n.shape[0], 1))
            translate_ = np.ones((n.shape[0], 1)) @ width_
            trans0 = np.hstack((translate_[:, 0, np.newaxis], zeros, zeros))
            trans1 = np.hstack((zeros, translate_[:, 1, np.newaxis], zeros))
            trans2 = np.hstack((zeros, zeros, translate_[:, 2, np.newaxis]))
            nodes_surr = np.vstack((nodes_surr, fipped_n[0] + trans0))  # add them
            nodes_surr = np.vstack((nodes_surr, fipped_n[0] - trans0))
            nodes_surr = np.vstack((nodes_surr, fipped_n[1] + trans1))
            nodes_surr = np.vstack((nodes_surr, fipped_n[1] - trans1))
            nodes_surr = np.vstack((nodes_surr, fipped_n[2] + trans2))
            nodes_surr = np.vstack((nodes_surr, fipped_n[2] - trans2))
            return nodes_surr, ni
        else:  # no nodes within
            return np.ones((0,)), np.ones((0,))

    def get_density(self, type:str, volumes:list = []):
        """ retrieves length, surface or volume density [cm/cm3, cm2/cm3, or cm3/cm3] per soil cells
            
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
        if not volumes:
            volumes = self.get_default_volumes_()
        ms = self.ms  # rename
        cell2seg = ms.cell2seg
        l = ms.segLength()
        radii = ms.getEffectiveRadii();
        n = int(ms.resolution.x * ms.resolution.y * ms.resolution.z)
        d = np.zeros((n,))
        for i in range(0, n):
            if i in cell2seg:
                for si in cell2seg[i]:
                    d[i] += f(radii[si], l[si])
        return np.divide(d, volumes)  # cm/cm3, cm2/cm3, or cm3/cm3

    def aggregate(self, seg_param):
        """ sums up params over the grid cells (parrams is given per segment)
        """
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
        """ averages params over the grid cells (parrams is given per segment)
        """
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

    def get_outer_radii(self, type:str, volumes:list = []):
        """ retrieves the outer radii of the perirhizal zones [cm] 
            
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
        else:
            print("PerirhizalPython.get_outer_radii() unknown type (should be 'length', 'surface', or 'volume')", type)
            raise

    def get_outer_radii_(self, type:str, volumes:list = None):
        """ python version of segOuterRadii, see get_outer_radii() """
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
        radii = ms.getEffectiveRadii();
        lengths = ms.segLength()
        width = ms.maxBound.minus(ms.minBound)
        outer_r = np.zeros((len(radii),))
        if not volumes:
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

    def get_outer_radii_bounded_voronoi(self):
        """ retrives the outer radii  [cm] and corresponding perirhizal volumes [cm3] 
            based on a 3D Voronoi diagram that are bounded by the soil volumes 
            using the approach of mirrored nodes. 
        """
        ms = self.ms  # rename

        min_b = ms.minBound
        max_b = ms.maxBound
        nodes_ = ms.nodes
        nodes = np.array([[n.x, n.y, n.z] for n in nodes_])  # to numpy array

        x = int(ms.resolution.x)
        y = int(ms.resolution.y)
        z = int(ms.resolution.z)

        width = np.array([max_b.x, max_b.y, max_b.z]) - np.array([min_b.x, min_b.y, min_b.z])
        print("making periodic", width)
        print("nodes", nodes.shape)
        print("resolution", x, y, z)
        nodes = self.make_periodic_(nodes, width)

        vol = np.empty((nodes.shape[0]))
        vol[:] = np.nan

        for i in range(0, x):
            for j in range(0, y):
                for k in range(0, z):

                    # print("\nIteration", i, j, k)

                    n, ni = self.mirror(nodes, i, j, k)
                    nn = ni.shape[0]

                    # print("ni", ni.shape, "vol", vol.shape)
                    # print("number of mirrored nodes", n.shape)
                    # print("n", n)
                    # print("ni", ni.shape)
                    # print(ni)
                    # dd

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
        radii = self.ms.getEffectiveRadii();
        lengths = self.ms.segLength()
        for i, v in enumerate(vol[1:]):  # seg_index = node_index -1
            if v > 0:
                outer_r[i] = np.sqrt(v / (np.pi * lengths[i]) + radii[i] * radii[i])

        return outer_r

    def get_outer_radii_voronoi(self):
        """ retrives the outer radii  [cm] and corresponding perirhizal volumes [cm3] 
            based on a 3D Voronoi diagram
        """
        min_b = self.ms.minBound
        max_b = self.ms.maxBound
        domain = pb.SDF_Cuboid(min_b, max_b)
        nodes_ = self.ms.nodes
        nodes = np.array([[n.x, n.y, n.z] for n in nodes_])  # to numpy array
        width = np.array([max_b.x, max_b.y, max_b.z]) - np.array([min_b.x, min_b.y, min_b.z])
        print("making periodic", width)
        nodes = self.make_periodic_(nodes, width)
        vor = Voronoi(nodes)
        vol = np.zeros(vor.npoints)

        for reg_num in vor.point_region:
            indices = vor.regions[reg_num]
            i = reg_num - 1
            if i >= 0:
                if -1 in indices:  # some regions can be opened
                    vol[i] = np.nan
                else:
                    inside = True  # convex polygon is within bounding box
                    for j in indices:
                        if domain.dist(vor.vertices[j]) >= 0:  # negative values mean inside the domain (sdf)
                            inside = False
                            break

                    if inside:
                        vol[i] = ConvexHull(vor.vertices[indices]).volume
                    else:
                        vol[i] = np.nan

        outer_r = vol.copy()
        radii = self.ms.getEffectiveRadii();
        lengths = self.ms.segLength()
        for i, v in enumerate(vol[1:]):  # seg_index = node_index -1
            if v > 0:
                outer_r[i] = np.sqrt(v / (np.pi * lengths[i]) + radii[i] * radii[i])

        return outer_r

    def get_voronoi_mesh(self, domain = None):

        nodes_ = self.ms.nodes  # make nodes periodic and pass all
        nodes = np.array([[n.x, n.y, n.z] for n in nodes_])  # to numpy array
        min_b = self.ms.minBound
        max_b = self.ms.maxBound
        width = np.array([max_b.x, max_b.y, max_b.z]) - np.array([min_b.x, min_b.y, min_b.z])
        nodes = self.make_periodic_(nodes, width)

        self = get_voronoi_mesh_(nodes, len(nodes))

    def get_node_mesh(self, nodes):
        grid = vtk.vtkUnstructuredGrid()
        points_array = vtk.vtkPoints()
        for n in nodes:
            points_array.InsertNextPoint(n)
        grid.SetPoints(points_array)
        return grid

    def get_voronoi_mesh_(self, vor_nodes, noc, crop_domain = None):
        """
        Creates a VTK grid containing all nodes, and only cells inside of domain 
        """
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

    def make_periodic_(self, nodes, width):
        """ maps the point periodically into [-width / 2., width / 2.] for x and y         
        """
        nodes_ = nodes.copy()
        for n in nodes_:
            n[0] = n[0] - ((n[0] + width[0] / 2.) // width[0]) * width[0]
            n[1] = n[1] - ((n[1] + width[1] / 2.) // width[1]) * width[1]
            # if n[0] < -width[0] / 2:
            #     print("n[0] < -width[0] / 2", n)
            # if n[0] > width[0] / 2:
            #     print("n[0] > width[0] / 2", n)
            # if n[1] < -width[1] / 2:
            #     print("n[1] < -width[0] / 2", n)
            # if n[1] > width[1] / 2:
            #     print("n[1] > width[0] / 2", n)
        return nodes_

    def to_range_(self, x:list, other:list = [], min_:float = -np.inf, max_:float = np.inf):
        """ returns ths list with x values within min_ and max_ (and drops nans) 
            drops the same entries from @param other list
        """
        y, z = [], []
        for i, x_ in enumerate(x):
            if (not np.isnan(x_)) and x_ >= min_ and x_ <= max_:
                y.append(x_)
                if other:
                    z.append(other[i])
        return y, z

    def get_default_volumes_(self):
        """ returns the soil volumes of a rectangular grid (default) """
        ms = self.ms  # rename
        width = ms.maxBound.minus(ms.minBound)
        vol = width.x * width.y * width.z / ms.resolution.x / ms.resolution.y / ms.resolution.z
        return np.ones((int(ms.resolution.x * ms.resolution.y * ms.resolution.z),)) * vol

    def flip_(self, nodes, center, axis):
        """ flips the nodes around the center according axis """
        n_ = nodes - np.ones((nodes.shape[0], 1)) @ center
        n_[:, axis] = -n_[:, axis]
        n_ = n_ + np.ones((nodes.shape[0], 1)) @ center
        return n_


if __name__ == "__main__":

    sand = [0.045, 0.43, 0.15, 3, 1000]
    loam = [0.08, 0.43, 0.04, 1.6, 50]
    clay = [0.1, 0.4, 0.01, 1.1, 10]

    hydrus_loam = [0.078, 0.43, 0.036, 1.56, 24.96]
    hydrus_clay = [0.068, 0.38, 0.008, 1.09, 4.8]
    hydrus_sand = [0.045, 0.43, 0.145, 2.68, 712.8]
    hydrus_sandyloam = [0.065, 0.41, 0.075, 1.89, 106.1]

    filename = "hydrus_loam"
    # sp = vg.Parameters(hydrus_loam)
    # vg.create_mfp_lookup(sp)
    peri = PerirhizalPython()
    # peri.create_lookup(filename, sp)  # takes some hours
    # peri.open_lookup(filename)

    peri.set_soil(vg.Parameters(loam))
    a = 0.1  # cm
    kr = 1.73e-4  # [1/day]
    rx = -15000  # cm
    sx = 0.  # cm
    rho = 1 / a
    inner_kr = a * kr
    rsx = peri.soil_root_interface_potentials([rx], [sx], [inner_kr], [rho])
    print("root soil interface", rsx, "cm")
    print("results into a flux of", kr * 2 * a * np.pi * (rsx - rx), "cm3/day")

