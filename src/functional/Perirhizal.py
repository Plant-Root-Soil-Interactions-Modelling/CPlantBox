import sys; sys.path.append(".."); sys.path.append("../..");

import plantbox as pb
from plantbox import Perirhizal
from plantbox import MappedSegments

import numpy as np


class PerirhizalPython(Perirhizal):
    """
    Wraps MappedSegments (or specialisations MappedPlant, MappedRootsystem)
    and adds functions to retrieve information on the perirhizal zones of single segments.
    """

    def __init__(self, ms:MappedSegments):
        """ @param ms is wrapped by this class """
        self.ms = ms

    def get_default_volumes_(self):
        """ returns the volumes of a rectangular grid (default) """
        width = ms.maxBound.minus(ms.minBound)
        vol = width.x * width.y * width.z / ms.resolution.x / ms.resolution.y / ms.resolution.z
        return np.ones((ms.resolution.x * ms.resolution.y * ms.resolution.z,)) * vol

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
        cell2seg = self.cell2seg
        l = self.segLength()
        n = np.max(cell2seg.keys())
        d = np.zeros((n,))
        for i in range(0, n):
            if i in cell2seg:
                for si in cell2seg[i]:
                    d[i] += f(self.radii[si], l[si])
        return np.divide(d, volumes)  # cm/cm3, cm2/cm3, or cm3/cm3

    def get_outer_radii(self, type:str, volumes:list = []):
        """ todo """
        if type == "length":
            return super().segOuterRadii(2, volumes)
        elif type == "surface":
            return super().segOuterRadii(1, volumes)
        elif type == "volume":
            return super().segOuterRadii(0, volumes)
        else:
            print("PerirhizalPython.get_outer_radii() unknown type (should be 'length', 'surface', or 'volume')", type)
            raise

    def get_outer_radii_voronoi(self, type):
        """ todo """
        raise

