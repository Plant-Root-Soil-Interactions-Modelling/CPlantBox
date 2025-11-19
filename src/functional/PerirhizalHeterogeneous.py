import sys; sys.path.append(".."); sys.path.append("../..");

import plantbox as pb
from Perirhizal import PerirhizalPython

import plantbox.functional.van_genuchten as vg

import numpy as np


class PerirhizalHetereogeneous(PerirhizalPython):
    """ 
    Same as PerirhizalPython but with multiple VG sets for heterogeneous soil layers (1D layers only)
    
    similar to cpp/soil_richards/richardsparams.hh but only for 1D
    
    TODO not finished yet 
    """

    def __init__(self, ms = None):
        """  ms      reference to MappedSegments """
        super().__init__(ms)

        self.lookup_tables = []  # optional array of 4d look up table to find soil root interface potentials
        self.sp = []  # corresponding array of van Genuchten soil parameter

        self.corners = [0]  # layers are given by there edges, i.e. [0, -10, -20] are two layers [0,-10] and [-10, -20]
        self.layer_indices = []  # van Genuchten set index per layer

        self.ind_ = []  # layer index per segment mid z coordinate

    def set_soil(self, sp):
        raise "PerirhizalHetereogeneous: use set_soils "

    def open_lookup(self, filename):
        raise "PerirhizalHetereogeneous: use open_lookup_tabless "

    def addLayers(self, corners, layer_indices):
        """ adds the layer geometry,  
            where soil parameter set with index @param layer_index applies, see addVanGenuchtenParamters, or addLookUpTable
        """
        I = np.argsort(corners)
        I = I[::-1]  # descending order
        scorners = np.array(corners)[I]
        slayer_indices = np.array(layer_indices)[I]
        self.corners = scorners
        self.layer_indices = slayer_indices

    def addVanGenuchtenParamters(self, sp_):
        """ 
        """
        self.sp.append(sp_)
        self.lookup_tables.append(None)

    def addLookUpTable(self, table_name):
        """
        """
        npzfile = np.load(table_name + ".npz")
        interface = npzfile["interface"]
        rx_, sx_, akr_, rho_ = npzfile["rx_"], npzfile["sx_"], npzfile["akr_"], npzfile["rho_"]
        soil = npzfile["soil"]
        table = RegularGridInterpolator((rx_, sx_, akr_, rho_), interface)  # method = "nearest" fill_value = None , bounds_error=False
        self.sp.append(vg.Parameters(soil))
        self.lookup_tables.append(table)

    def test(self):
        """ checks if sizes are correct """
        assert len(self.boxes_min) == len(self.boxes_max) == len(self.layer_indices), "PerirhizalHetereogeneous.test(): Lengths should agree"
        assert len(self.sp_) == len(self.lookup_tables), "PerirhizalHetereogeneous.test(): Lengths should agree"
        assert len(self.lookup_tables) > np.max(self.layer_indices), "PerirhizalHetereogeneous.test(): there is a layer index exceeding the table"

    def create_layer_indices(self):
        """
        """
        print("create_layer_indices")
        z = np.array(self.ms.getSegmentZ())
        self.ind_ = np.zeros(z.shape, dtype = int)
        for i in range(1, len(self.corners)):  # per layer
            self.ind_ += (np.logical_and(z <= np.array(self.corners[i - 1]), z > np.array(self.corners[i]))) * (i - 1)
        self.ind_ = np.array(list(self.layer_indices[i] for i in self.ind_))

    def soil_root_interface_potentials(self, rx, sx, inner_kr, rho):
        """
        finds matric potentials at the soil root interface for as all segments
        uses a look up tables if present (see create_lookup, and open_lookup) 
        addVanGenuchtenDomain
        rx             xylem matric potential [cm]
        sx             bulk soil matric potential [cm]
        inner_kr       root radius times hydraulic conductivity [cm/day] 
        rho            geometry factor [1] (outer_radius / inner_radius)
        """
        assert len(rx) == len(sx) == len(inner_kr) == len(rho), "rx, sx, inner_kr, and rho must have the same length"
        if not len(ind_) == len(x):
          self. create_layer_indices(self)
        assert len(rx) == len(ind_), "number of segments and xylem potentials must have the same length"

        rsx = np.zeros(rx.shape)

        for i in range(0, np.max(self.layer_indices)):  # per soil parameter set
            self.lookup_table = self.lookup_tables[i]  # used for table look up
            if self.lookup_tables:  # if it was set
                i0 = self.ind_ == i  # binary indexing of the segments with layer i
                rsx[i0] = self.soil_root_interface_potentials_table(rx[i0], sx[i0], inner_kr[i0], rho[i0])
            else:
                rsx[i0] = np.array([PerirhizalPython.soil_root_interface_(rx[i], sx[i], inner_kr[i], rho[i], self.sp) for i in range(0, len(rx))])
                rsx = rsx[:, 0]

        return rsx


if __name__ == "__main__":

    plant = pb.MappedPlant()
    path = "../../modelparameter/structural/rootsystem/"
    name = "Anagallis_femina_Leitner_2010"
    plant.readParameters(path + name + ".xml")
    plant.initialize()
    plant.simulate(4, False)
    # setGrid

    peri = PerirhizalHetereogeneous(plant)

    corners = [0, -2, -4, -8, -120, -150, -180, -190]  # 7 layers
    peri.addLayers(corners, range(0, 8))

    peri.create_layer_indices()
    print(peri.ind_)

    print("fin")
