import sys; sys.path.append(".."); sys.path.append("../..");

import plantbox as pb
from plantbox import Perirhizal
from plantbox import MappedSegments

import numpy as np
from scipy.spatial import ConvexHull, Voronoi
import scipy.linalg as la
import vtk


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
        min_ = np.array([ms.minBound.x + i * cell_width[0], ms.minBound.y + j * cell_width[1], ms.minBound.z + k * cell_width[2]])
        max_ = min_ + cell_width
        return min_, max_

    def nodes_within(self, nodes, min_, max_):
        """ returns the nodes inside the bounding box given by min_, and max_ """
        one = np.ones((nodes.shape[0], 1))
        n_, ni = [], []
        for i, n in enumerate(nodes):
            if n[0] >= min_[0] and n[1] >= min_[1] and n[2] >= min_[2] and n[0] < max_[0] and n[1] < max_[1] and n[2] < max_[2]:
                n_.append(n)
                ni.append(i)
        # print("nodes within", min_, max_, ":", len(n_))
        return np.array(n_), np.array(ni)

    def flip_(self, nodes, center, axis):
        """ flips the nodes around the center according axis """
        n_ = nodes - np.ones((nodes.shape[0], 1)) @ center
        n_[:, axis] = -n_[:, axis]
        n_ = n_ + np.ones((nodes.shape[0], 1)) @ center
        return n_

    def mirror(self, nodes, i:int, j:int, k:int):
        """ adds mirrored nodes to the 6 sides of the cube with index i, j, k"""
        min_, max_ = self.get_bounds(i, j, k)
        width_ = max_ - min_
        width_ = np.expand_dims(width_, axis = 0)
        center_ = min_ + width_ / 2.
        n, ni = self.nodes_within(nodes, min_, max_)
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

    def aggregate(self, seg_param):
        """ sums up params over the grid cells (parrams is given per segment)
        
        TODO move to MappedSegments
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
        
        TODO move to MappedSegments
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
                                    print("PerirhizalPython.get_outer_radii_bounded_voronoi(): nan encountered")
                                else:
                                    vol[ni[idx]] = ConvexHull(vor.vertices[indices]).volume
                                    # print(ni[i_])
                                    # print("index", ni[i_], "vol", vol[ni[i_]])
                            # else:
                            #     if i_ < 0:
                            #         print("When does that happen, again?", idx, nn, indices)

        if np.sum(np.isnan(vol)) > 1:
            print("PerirhizalPython.get_outer_radii_bounded_voronoi(): nan encountered", np.sum(np.isnan(vol)), np.where(np.isnan(vol)))

        outer_r = np.zeros((vol.shape[0] - 1,))
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

