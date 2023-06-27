import sys; sys.path.append(".."); sys.path.append("../..");

import plantbox as pb
from plantbox import Perirhizal
from plantbox import MappedSegments

import numpy as np
from scipy.spatial import ConvexHull, Voronoi
import scipy.linalg as la
import vtk
from mercurial.templatefuncs import min_, max_


class PerirhizalPython(Perirhizal):
    """
    Retrieve information on the perirhizal zones
    
    Wraps MappedSegments (or specialisations MappedPlant, MappedRootsystem)
    and adds functions to retrieve information on the perirhizal zones of single segments.
    """

    def __init__(self, ms:MappedSegments):
        """ @param ms is wrapped by this class """
        super().__init__(ms)

    def get_default_volumes_(self):
        """ returns the soil volumes of a rectangular grid (default) """
        ms = self.ms  # rename
        width = ms.maxBound.minus(ms.minBound)
        vol = width.x * width.y * width.z / ms.resolution.x / ms.resolution.y / ms.resolution.z
        return np.ones((int(ms.resolution.x * ms.resolution.y * ms.resolution.z),)) * vol

    def get_bounds(self, i:int, j:int, k:int):
        """ returns the bounding box (minx ,miny ,minz, maxx, maxy, maxz) of soil cell with index i, j, k of a rectangular grid """
        ms = self.ms  # rename
        width = ms.maxBound.minus(ms.minBound)
        cell_width = np.array([width.x / ms.resolution.x, width.y / ms.resolution.y, width.z / ms.resolution.z])
        min_ = i * cell_width
        max_ = (i + 1) * cell_width
        return min_, max_

    def nodes_within(self, min_, max_):
        """ returns the nodes inside the bounding box given by min_, and max_ """
        nodes_ = self.ms.nodes
        nodes = np.array([[n.x, n.y, n.z] for n in nodes_])
        one = np.ones((nodes.shape[0], 1))
        n_, ni = [], []
        for i, n in enumerate(nodes):
            if n[0] >= min_[0] and n[1] >= min_[1] and n[2] >= min_[2] and n[0] < max_[0] and n[1] < max_[1] and n[2] < max_[2]:
                n_.append(n)
                ni.append(i)
        print("nodes within", min_, max_, ":", len(n_))
        print(n_)
        return np.array(n_), np.array(ni)

    def flip_(self, nodes, center, axis):
        """ flips the nodes around the center according axis """
        n_ = nodes - np.ones((nodes.shape[0], 1)) @ center
        n_[:, axis] = -n_[:, axis]
        n_ = n_ + center
        return n_

    def mirror_(self, i:int, j:int, k:int):
        """ adds mirrored nodes to the 6 sides of the cubes"""
        min_, max_ = self.get_bounds(i, j, k)
        width_ = max_ - min_
        width_ = np.expand_dims(width_, axis = 0)
        center_ = min_ + width_ / 2.
        n, ni = self.nodes_within(min_, max_)
        if n.shape[0] > 0:
            nodes_surr = n
            fipped_n = [self.flip_(n, center_, i_) for i_ in range(0, 3)]  # flip them ...
            translate_ = np.ones((n.shape[0], 1)) @ width_
            print(width_.shape)
            print(n.shape)
            print("translate_", translate_.shape)
            for i_ in range(0, 3):
                nodes_surr = np.vstack((nodes_surr, fipped_n[i_] + translate_[:, i_]))  # add them
                nodes_surr = np.vstack((nodes_surr, fipped_n[i_] - translate_[:, i_]))
            return nodes_surr, ni
        else:  # no nodes within
            return np.ones((0,)), np.ones((0,))

    def get_density(self, type:str, volumes:list = []):
        """ retrives length, surface or volume density [cm/cm3, cm2/cm3, or cm3/cm3]
            
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
        radii = ms.radii
        n = int(ms.resolution.x * ms.resolution.y * ms.resolution.z)
        d = np.zeros((n,))
        for i in range(0, n):
            if i in cell2seg:
                for si in cell2seg[i]:
                    d[i] += f(radii[si], l[si])
        return np.divide(d, volumes)  # cm/cm3, cm2/cm3, or cm3/cm3

    def get_outer_radii(self, type:str, volumes:list = []):
        """ retrives the outer radii of the perirhizal zones [cm]
            
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

    def get_outer_radii_(self, type:str, volumes:list = []):
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
        radii = ms.radii
        lengths = ms.segLength()
        width = ms.maxBound.minus(ms.minBound)
        outer_r = np.zeros((len(radii),))
        if not volumes:
            n = ms.resolution.x * ms.resolution.y * ms.resolution.z
            volumes = np.ones(int(n)) * (width.x * width.y * width.z) / n
            print("PerirhizalPython.get_outer_radii_(): each soil volume has", (width.x * width.y * width.z) / n, "cm3")

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
        width = np.array([max_b.x, max_b.y, max_b.z]) - np.array([min_b.x, min_b.y, min_b.z])
        print("making periodic", width)
        print("nodes", nodes.shape)

        nodes = self.make_periodic_(nodes, width)

        x = int(ms.resolution.x)
        y = int(ms.resolution.y)
        z = int(ms.resolution.z)

        vol = np.zeros((nodes.shape[0]))

        for i in range(0, x):
            for j in range(0, y):
                for k in range(0, z):
                    n, ni = self.mirror_(i, j, k)
                    nn = ni.shape[0]
                    # print("n", n)
                    # print("ni", ni)
                    if nn > 0:
                        vor = Voronoi(n)

                        for reg_num in vor.point_region:
                            indices = vor.regions[reg_num]
                            i_ = reg_num - 1
                            if i_ >= 0 and i_ < nn:
                                if -1 in indices:  # some regions can be opened
                                    print(i_)
                                    print(ni)
                                    vol[ni[i_]] = np.nan
                                else:
                                    vol[ni[i_]] = ConvexHull(vor.vertices[indices]).volume
                            else:
                                if i_ < 0:
                                    print("When does that happen, again?")

        outer_r = vol.copy()
        radii = self.ms.radii
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
        radii = self.ms.radii
        lengths = self.ms.segLength()
        for i, v in enumerate(vol[1:]):  # seg_index = node_index -1
            if v > 0:
                outer_r[i] = np.sqrt(v / (np.pi * lengths[i]) + radii[i] * radii[i])

        return outer_r

    def get_voronoi_mesh(self, domain = None):
        """
        Creates a VTK grid containing all nodes, and only cells inside of domain 
        """
        min_b = self.ms.minBound
        max_b = self.ms.maxBound
        width = np.array([max_b.x, max_b.y, max_b.z]) - np.array([min_b.x, min_b.y, min_b.z])
        if not domain:
            domain = pb.SDF_Cuboid(min_b, max_b)

        nodes_ = self.ms.nodes
        nodes = np.array([[n.x, n.y, n.z] for n in nodes_])  # to numpy array
        print(width)
        print("*")
        nodes = self.make_periodic_(nodes, width)
        vor = Voronoi(nodes)

        grid = vtk.vtkUnstructuredGrid()  # Create a VTK unstructured grid

        points_array = vtk.vtkPoints()  # Add the Voronoi vertices as points
        for v in vor.vertices:
            if domain.dist(v) <= 0:
                points_array.InsertNextPoint(v)
            else:
                points_array.InsertNextPoint([0, 0, 0])  # simpler for vizualizing in ParaView
        grid.SetPoints(points_array)

        vol = []
        for reg_num in vor.point_region:
            indices = vor.regions[reg_num]
            i = reg_num - 1
            if i >= 0:
                if not -1 in indices:  # not an open region
                    inside = True
                    for j in indices:
                        if domain.dist(vor.vertices[j]) >= 0:  # negative values mean inside the domain (sdf)
                            inside = False
                            break
                    if inside:
                        vol.append(ConvexHull(vor.vertices[indices]).volume)
                        id_array = vtk.vtkIdList()
                        for vertex_index in indices:
                            id_array.InsertNextId(vertex_index)
                        grid.InsertNextCell(vtk.VTK_CONVEX_POINT_SET, id_array)

        cell_id = vtk.vtkDoubleArray()
        cell_id.SetName("cell_id")
        cell_id.SetNumberOfValues(len(vol))
        for j in range(0, len(vol)):
            cell_id.SetValue(j, j)
        celldata = grid.GetCellData()
        celldata.AddArray(cell_id)

        vols = vtk.vtkDoubleArray()
        vols.SetName("volumes")
        vols.SetNumberOfValues(len(vol))
        for j in range(0, len(vol)):
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

    def to_range_(self, x:list, min_:float, max_:float):
        """ returns ths list with values within min_ and max_ (and drops nans) """
        y = []
        for x_ in x:
            if (not np.isnan(x_)) and x_ >= min_ and x_ <= max_:
                y.append(x_)
        return y

