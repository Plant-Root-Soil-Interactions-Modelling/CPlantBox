import sys
import numpy as np

# if you need path tweaks, you can keep these two, but NOT an import of this file itself
# sys.path.append("..")
# sys.path.append("../..")

import plantbox as pb

from functional.Perirhizal import PerirhizalPython
import functional.van_genuchten as vg
from scipy.interpolate import RegularGridInterpolator

# global cache so multiple plants can reuse the same table
_TABLE_CACHE = {}  # key: filename (str) -> (table, vg.Parameters)



class PerirhizalHeterogeneous(PerirhizalPython):
    """
    Perirhizal with several soil/lookup sets, chosen by z-layers.
    """

    def __init__(self, ms=None):
        super().__init__(ms)

        # parallel lists: same index = same soil
        self.lookup_tables = []   # list of RegularGridInterpolator or None
        self.sp_list = []         # list of vg.Parameters

        # layer description
        self.corners = [0.0]      # z values, e.g. [0, -25, -50]
        self.layer_indices = []   # which soil index for each interval

        self.boxes = []
        self.seg_xyz = None   # <- optional, only set if user gives it

    def add_box(self, p1, p2, soil_idx):
        self.boxes.append((np.array(p1), np.array(p2), soil_idx))

    def set_segment_xyz(self, xyz):
        """xyz: shape (nseg, 3), segment midpoints in soil coords"""
        self.seg_xyz = np.asarray(xyz)
        
    def add_vg_parameters(self, sp_):
        # register for fast_mfp so PerirhizalPython.soil_root_interface_ can use it
        vg.create_mfp_lookup(sp_)
        self.sp_list.append(sp_)
        self.lookup_tables.append(None)

    def add_lookup_table(self, filename):
        """
        Load a precomputed .npz lookup table and register it as a soil.
        Uses a module-level cache so the same file is only loaded once,
        even if multiple plants / perirhizal objects call this.
        """
        global _TABLE_CACHE

        # 1) seen before? just reuse
        if filename in _TABLE_CACHE:
            table, sp = _TABLE_CACHE[filename]
            self.lookup_tables.append(table)
            self.sp_list.append(sp)
            return

        # 2) otherwise load from disk
        npzfile = np.load(filename + ".npz")
        interface = npzfile["interface"]
        rx_, sx_, akr_, rho_ = (
            npzfile["rx_"],
            npzfile["sx_"],
            npzfile["akr_"],
            npzfile["rho_"],
        )
        soil = npzfile["soil"]

        table = RegularGridInterpolator((rx_, sx_, akr_, rho_), interface)

        # build vg.Parameters and register its mfp
        sp = vg.Parameters(soil)
        vg.create_mfp_lookup(sp)

        # put into global cache
        _TABLE_CACHE[filename] = (table, sp)

        # and also into THIS instance
        self.lookup_tables.append(table)
        self.sp_list.append(sp)


    # ------------------------------------------------------
    # layers
    # ------------------------------------------------------
    def add_layers(self, corners, layer_indices):
        """
        corners: e.g. [0, -25, -50]
        layer_indices: e.g. [0, 1]  (same length as corners-1)
        """
        # sort descending in z
        I = np.argsort(corners)[::-1]
        corners = np.array(corners)[I]
        layer_indices = np.array(layer_indices)[I[:-1]]

        self.corners = corners
        self.layer_indices = layer_indices

    def create_layer_indices(self):
        # always do z-layers first
        z = np.array(self.ms.getSegmentZ())
        seg_soil = np.zeros_like(z, dtype=int)

        for i in range(len(self.corners) - 1):
            upper = self.corners[i]
            lower = self.corners[i + 1]
            mask = (z <= upper) & (z > lower)
            seg_soil[mask] = self.layer_indices[i]

        # if we got full xyz, we can cut out boxes and overwrite
        if self.seg_xyz is not None and len(self.boxes) > 0:
            xyz = self.seg_xyz
            for (p1, p2, soil_idx) in self.boxes:
                m = (
                    (xyz[:, 0] >= p1[0]) & (xyz[:, 0] <= p2[0]) &
                    (xyz[:, 1] >= p1[1]) & (xyz[:, 1] <= p2[1]) &
                    (xyz[:, 2] >= p1[2]) & (xyz[:, 2] <= p2[2])
                )
                seg_soil[m] = soil_idx

        self.seg_soil_index = seg_soil

    def set_segment_soil_index(self, seg_soil_index: np.ndarray):
        self.seg_soil_index = np.asarray(seg_soil_index, dtype=int)


    # ------------------------------------------------------
    # main interface
    # ------------------------------------------------------
    def soil_root_interface_potentials(self, rx, sx, inner_kr, rho):
        """
        Same signature as base class but now chooses soil/table per segment.
        We convert inputs to numpy arrays so boolean masks work.
        """
        # make sure we can mask
        rx = np.asarray(rx)
        sx = np.asarray(sx)
        inner_kr = np.asarray(inner_kr)
        rho = np.asarray(rho)

        # (re)build segment → soil mapping if needed
        if self.seg_soil_index is None or len(self.seg_soil_index) != len(rx):
            self.create_layer_indices()

        rsx = np.zeros_like(rx, dtype=float)

        # loop over all soil sets we have stored
        for soil_idx in range(len(self.sp_list)):
            mask = (self.seg_soil_index == soil_idx)
            if not np.any(mask):
                continue

            table = self.lookup_tables[soil_idx]
            if table is not None:
                # use lookup-table path
                self.lookup_table = table  # reuse parent’s helper
                rsx[mask] = self.soil_root_interface_potentials_table(
                    rx[mask], sx[mask], inner_kr[mask], rho[mask]
                )
            else:
                # fall back to analytical VG for this soil
                sp = self.sp_list[soil_idx]
                idxs = np.where(mask)[0]
                vals = []
                for j in idxs:
                    vals.append(
                        PerirhizalPython.soil_root_interface_(
                            rx[j], sx[j], inner_kr[j], rho[j], sp
                        )
                    )
                rsx[mask] = np.array(vals)

        return rsx

