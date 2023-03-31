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

    def get_outer_radii_voronoi(self):
        """ retrives the outer radii  [cm] and corresponding perirhizal volumes [cm3] 
            based on a 3D Voronoi diagram
        """
        min_b = self.ms.minBound
        max_b = self.ms.maxBound
        domain = pb.SDF_Cuboid(min_b, max_b)
        nodes_ = self.ms.nodes
        nodes = np.array([[n.x, n.y, n.z] for n in nodes_])  # to numpy array
        vor = Voronoi(nodes)
        vol = np.zeros(vor.npoints)

        for reg_num in vor.point_region:
            indices = vor.regions[reg_num]
            i = reg_num - 1
            if i >= 0:
                if -1 in indices:  # some regions can be opened
                    vol[i] = np.nan
                else:
                    inside = True
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

    def intersect_polygon_quader_(self, polygon, quader):
        """ 
            This code first defines the polygon and quader vertices as numpy arrays. 
            Then, it defines the planes of the polygon and the quader using the cross product and dot product. 
            It then checks for an intersection using the separating axis theorem by computing the dot product of each vertex with each plane and checking
        """
        # Define planes of polygon
        planes_polygon = []
        for i in range(len(polygon)):
            p1 = polygon[i]
            p2 = polygon[(i + 1) % len(polygon)]
            normal = np.cross(p2 - p1, [0, 0, 1])
            d = -np.dot(normal, p1)
            planes_polygon.append(np.hstaverck([normal, d]))

        # Check for intersection using separating axis theorem
        for plane in planes_polygon + planes_quader:
            dpoly = np.dot(polygon, plane[:3]) + plane[3]
            dquader = np.dot(quader, plane[:3]) + plane[3]
            if np.all(dpoly < 0) or np.all(dquader < 0):
                print("No intersection")
                break
        else:
            # Find intersection points
            intersection_points = []
            for i in range(len(polygon)):
                p1 = polygon[i]
                p2 = polygon[(i + 1) % len(polygon)]
                for plane in planes_quader:
                    normal = plane[:3]
                    d = plane[3]
                    t = -(np.dot(normal, p1) + d) / np.dot(normal, p2 - p1)
                    if 0 <= t <= 1:
                        intersection_points.append(p1 + t * (p2 - p1))
            return intersection_points

    def to_range_(self, x:list, min_:float, max_:float):
        """ returns ths list with values within min and max (and drops nans) """
        y = []
        for x_ in x:
            if (not np.isnan(x_)) and x_ >= min_ and x_ <= max_:
                y.append(x_)
        return y

