import sys; sys.path.append(".."); sys.path.append("../..");

import plantbox as pb
from plantbox import Perirhizal
from plantbox import MappedSegments

import numpy as np
from scipy.spatial import ConvexHull, Voronoi
import scipy.linalg as la


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
        """ returns the volumes of a rectangular grid (default) """
        ms = self.ms  # rename
        width = ms.maxBound.minus(ms.minBound)
        vol = width.x * width.y * width.z / ms.resolution.x / ms.resolution.y / ms.resolution.z
        return np.ones((int(ms.resolution.x * ms.resolution.y * ms.resolution.z),)) * vol

    def get_density(self, type:str, volumes:list = []):
        """ retrives length, surface or volume density
            volumes of the finite volume cells, 
            if list is empty a recangular grid is assumed (its data is stored in MappedSegments), 
            see get_default_volumes_() """
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
        """ todo """
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
        """ python version of segOuterRadii """
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
            # print("PerirhizalPython.get_outer_radii_() each volume has", (width.x * width.y * width.z) / n, "cm3")

        for cell_id, seg_ids in cell2seg.items():
            tt = np.sum(np.array([f(radii[i], lengths[i]) for i in seg_ids]))
            for si in seg_ids:
                t = f(radii[si], lengths[si]) / tt  # proportionality factor (must sum up to == 1 over cell)
                v = t * volumes[cell_id]  # target volume
                outer_r[si] = np.sqrt(v / (np.pi * lengths[si]) + radii[si] * radii[si])

        return outer_r

    def get_outer_radii_voronoi(self):
        """ todo """

        min_b = self.ms.minBound
        max_b = self.ms.maxBound
        domain = pb.SDF_Cuboid(min_b, max_b)

        nodes_ = self.ms.nodes
        nodes = np.array([[n.x, n.y, n.z] for n in nodes_])
        segs = self.ms.segments
        lengths = self.ms.segLength()

        vor = Voronoi(nodes)
        vol = np.zeros(vor.npoints)
        outer_r = np.zeros(vol.shape)
        # quader, planes_quader = self.create_quader_from_bounding_box_([min_b.x, min_b.y, min_b.z], [max_b.x, max_b.y, max_b.z])
        # print(quader, planes_quader)

        for i, reg_num in enumerate(vor.point_region):
            # print(i, reg_num)
            indices = vor.regions[reg_num]
            if -1 in indices:  # some regions can be opened
                vol[i] = -1
            else:
                inside = True
                for j in indices:
                    if domain.dist(vor.vertices[j]) >= 0:  # negative values mean inside the domain (sdf)
                        inside = False
                        break
                if inside:
                    vol[i] = ConvexHull(vor.vertices[indices]).volume
                else:
                    vol[i] = 0

        print("shaping... ", vol.shape, nodes.shape, vor.npoints, len(vor.point_region), np.min(vor.point_region), np.max(vor.point_region))  # node indices
        # one point too much in point_region?????

        dd

        for si, seg in enumerate(segs):
            v = vol[seg.y - 1]
            outer_r[si] = np.sqrt(v / (np.pi * lengths[si]) + radii[si] * radii[si])

        # for i, reg_num in enumerate(vor.point_region):
        #     indices = vor.regions[reg_num]
        #     if -1 in indices:  # some regions can be opened
        #         vol[i] = 0
        #     else:
        #         vol[i] = ConvexHull(vor.vertices[indices]).volume
        #         vol[i] = min(vol[i], 8.)

        return outer_r  # TODO unfinished...

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
            planes_polygon.append(np.hstack([normal, d]))

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

    def get_voronoi_mesh(self):
        grid = vtkUnstructuredGrid()  # Create a VTK unstructured grid

        points_array = vtkPoints()  # Add the Voronoi vertices as points
        for v in ver:
            for j in range(0, 3):
                v[j] = min(v[j], max_b[j])
                v[j] = max(v[j], min_b[j])
            points_array.InsertNextPoint(v)
        grid.SetPoints(points_array)
        print("voronoi nodes: ", ver.shape)

        # print("vor.regions", vor.regions)
        cc = 0
        for region in vor.regions:
            if len(region) > 0 and (-1 not in region):
                cc += 1
                id_array = vtk.vtkIdList()
                id_array.InsertNextId(len(region))
                # print("added an id array with length", len(region))
                for vertex_index in region:
                    id_array.InsertNextId(vertex_index)
                grid.InsertNextCell(vtk.VTK_CONVEX_POINT_SET, id_array)

        cell_id = vtk.vtkDoubleArray()
        cell_id.SetName("cell_id")
        cell_id.SetNumberOfValues(cc)
        vc = 0
        for j in range(0, cc):
            cell_id.SetValue(vc, vc)
            vc += 1
        celldata = grid.GetCellData()
        celldata.AddArray(cell_id)

        cell_id = vtk.vtkDoubleArray()
        cell_id.SetName("volumes")
        cell_id.SetNumberOfValues(cc)
        vc = 0
        for j, region in enumerate(vor.regions):
            if len(region) > 0 and (-1 not in region):
                cell_id.SetValue(vc, vol0[j])
                vc += 1
        celldata = grid.GetCellData()
        celldata.AddArray(cell_id)

        # Write the VTK unstructured grid to a ParaView VTU file
        writer = vtkXMLUnstructuredGridWriter()
        writer.SetFileName('rootsystem.vtu')
        writer.SetInputData(grid)
        writer.Write()

