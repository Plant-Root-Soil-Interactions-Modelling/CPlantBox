import rsml.rsml_reader as rsml_reader

import numpy as np


class RsmlData:
    """ 
    Stores RSML data
    
    use open_rsml to open a rsml file 
    * automatically converts units to cm and day 
    * if necessary converts 2d -> 3d 
    * automatically detects xml tags for radii, creation times, and subTypes, and converts to per node    
    
    used by the different GUIS in CPlantBox/gui
    """

    def __init__(self):
        self.fname = None  # file name (str)
        self.polylines = None  # from rsml (list of list of list)
        self.properties = None  # from rsml (dict of list) value per root
        self.functions = None  # from rsml (dict of list of list) value per node
        self.metadata = None  # from rsml (rsml_writer.Metadata)
        self.radii = None  # list of list (per node)
        self.cts = None  # creation times, list of list (per node)
        self.max_ct = 0.  # maximal creation time
        self.types = None  # subTypesm, list of list (per node)
        self.tagnames = None  # list of str (original tag names, None if not found)

        # add more where needed: key = cm *value (for length), and key = day * value (for time)
        self.scales_ = {"pixel":1, "px":1, "dots": 1,
                   "cm": 1, "mm": 0.1, "dm": 10., "m": 100,
                   "h": 1. / 24., "s": 1 / (24.*3600), "sec": 1 / (24.*3600), "d": 1., "day": 1, "days": 1}

    def exists(self):
        """ true if a rsml was set """
        return self.polylines is not None

    def set_rsml(self, polylines, properties, functions, metadata):
        """ setter for rsml fields """
        self.polylines = polylines
        self.properties = properties
        self.functions = functions
        self.metadata = metadata

    def set_selected(self, radii, cts, types, tagnames):
        """ setter for selected tags """
        self.radii = radii
        if np.isnan(cts[0]):  # nan indicates creation times not given in rsml
            self.cts = np.zeros((len(cts),))
        else:
            self.cts = cts
        self.max_ct = np.max(self.cts)
        self.types = types
        self.tagnames = tagnames

    def open_rsml(self, fname, shift_z = False):
        """ opens an rsml file into self.data, using rsml_reader (in CPlantBox/src/python_modules)              
        converts units to cm and day 
        if necessary converts 2d -> 3d, 
        """
        polylines, properties, functions, metadata = rsml_reader.read_rsml(fname)
        # if metadata.software == "archisimple":
        #     print("DataModel.open_rsml(): special rules for archisimple: switch -y and z")
        #     new_polylines = []
        #     for pl in polylines:
        #         pl_ = []
        #         for p in pl:
        #             p_ = [p[0], p[2], -p[1]]
        #             pl_.append(p_)
        #         new_polylines.append(pl_)
        #     polylines = new_polylines

        print("DataModel.open_rsml(): scale to cm", metadata.scale_to_cm)
        self.set_rsml(polylines, properties, functions, metadata)
        self.scale_polylines_()  # converts units
        self.check_polylines_2d_(shift_z)  # 2d -> 3d
        radii, cts, types, tagnames = rsml_reader.get_parameter(polylines, functions, properties)  # paramter per node
        self.set_selected(radii, cts, types, tagnames)
        self.scale_selected_()  # converts units of special fields radii, cts, types

    def scale_polylines_(self):
        """  scales nodes, see rsml_writer.Metadata, and self.scale_to_cm 
        """
        scale = self.metadata.scale_to_cm  # default length scales
        for i in range(0, len(self.polylines)):
            for j in range(0, len(self.polylines[i])):
                for k in range(0, 3):
                    self.polylines[i][j][k] *= scale

    def check_polylines_2d_(self, shift_z = False):
        """  converts 2d image coordinates to 3d coordinates
            shift_z determines if the roots system seed is shifted to -3 cm
        """
        nodes, segs = rsml_reader.get_segments(self.polylines, self.properties)  # fetch nodes and segments
        maxz = np.max(nodes[:, 2])
        minz = np.min(nodes[:, 2])
        if maxz >= 0 and minz >= 0:  # image coordinates in px often start in lower left corner
            print("DataModel.check_polylines_2d_(): assuming image coordinates, y-centered and z-flipped ")
            miny = np.min(nodes[:, 1])
            yy = np.max(nodes[:, 1]) - miny
            for pl in self.polylines:  # both (pl and node) are references
                for node in pl:
                    node[2] = -node[2]
                    node[1] = node[1] - miny - yy / 2
        if shift_z:
            print("DataModel.check_polylines_2d_(): root system seed is shifted to (x,y,-3) ")
            z = self.polylines[0][0][2]
            print("z", z)
            for pl in self.polylines:  # both (pl and node) are references
                for node in pl:
                    node[2] = node[2] - z + (-3)

    def scale_selected_(self):
        """  scales radius and creation times, see rsml_writer.Metadata, and self.scale_to_cm
        shifts types, in a way they start at zero
        """
        scale = self.metadata.scale_to_cm  # default length scales
        # radii
        extra_str = ""
        if self.tagnames[0]:
            if self.tagnames[0] in self.metadata.properties:
                r_scale = self.scales_[self.metadata.properties[self.tagnames[0]].unit]
            else:  # assume same scaling as polylines
                r_scale = scale
        else:  # assume same scaling as polylines
            r_scale = scale
        if self.metadata.software == "smartroot":
            r_scale = scale
            extra_str = " (smartroot)"
        if self.metadata.software == "archisimple":
            r_scale = scale
            extra_str = " (archisimple)"
        print("DataModel.scale_selected_():radius length scale" + extra_str, r_scale)
        for i in range (0, len(self.radii)):
            self.radii[i] *= r_scale
        # creation times
        cts_scale = 1.  # assume it is in days
        if self.tagnames[1]:
            if self.tagnames[1] in self.metadata.properties:
                cts_scale = self.scales_[self.metadata.properties[self.tagnames[1]].unit]
                print("DataModel.scale_rsml() temporal scale", cts_scale)
        for i in range (0, len(self.cts)):
            self.cts[i] *= cts_scale
        # types
        min_types = np.min(self.types)
        for i in range (0, len(self.types)):
            self.types[i] -= min_types
