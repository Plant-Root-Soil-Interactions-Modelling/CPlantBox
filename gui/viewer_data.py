import sys; sys.path.append("../src/python_modules/"); sys.path.append("../")

import rsml_reader
import plantbox as pb

import numpy as np


class DataModel:
    """ MVC """

    def __init__(self):
        self.fname = None  # file name
        self.polylines = None  # from rsml
        self.properties = None  # from rsml
        self.functions = None  # from rsml
        self.metadata = None  # from rsml
        self.radii, self.cts, self.types, self.tagnames = None, None, None, None  # selected from rsml
        self.analyser = None  # created by convert_to_analyser
        self.mapped_segments = None  # created by convert_to_analyser
        self.max_ct = 0.  # created by convert_to_analyser
        # add more where needed, key = cm *value (for length), and key = day * value (for time)
        self.scales_ = {"pixel":1, "px":1, "dots": 1,
                   "cm": 1, "mm": 0.1, "dm": 10., "m": 100,
                   "h": 1. / 24., "s": 1 / (24.*3600), "sec": 1 / (24.*3600), "d": 1., "day": 1, "days": 1}

    def exists(self):
        """ true if a rsml was set """
        return self.polylines is not None

    def set_rsml(self, polylines, properties, functions, metadata):
        """ setter for rmsl fields """
        self.polylines = polylines
        self.properties = properties
        self.functions = functions
        self.metadata = metadata

    def set_selected(self, radii, cts, types, tagnames):
        """ setter for selected tags """
        self.radii = radii
        self.cts = cts
        self.types = types
        self.tagnames = tagnames

    def open_rsml(self, fname):
        """ opens an rsml file into self.data """
        polylines, properties, functions, metadata = rsml_reader.read_rsml(fname)
        self.set_rsml(polylines, properties, functions, metadata)
        radii, cts, types, tagnames = rsml_reader.get_parameter(polylines, functions, properties)  # paramter per node
        self.set_selected(radii, cts, types, tagnames)
        print("DataModel.open_rsml(): scale to cm", metadata.scale_to_cm)
        self.scale_rsml()
        self.convert_to_analyser()

    def scale_rsml(self):
        """ 
        automatically scales coordinates and creation times, see rsml_writer.Metadata, and self.scales_
        """
        if self.exists():
            scale = self.metadata.scale_to_cm  # default length scales
            # radii
            if self.tagnames[0]:
                if self.tagnames[0] in self.metadata.properties:
                    r_scale = self.scales_[self.metadata.properties[self.tagnames[0]].unit]
                    print("radius length scale", r_scale)
                else:  # assume same scaling as polylines
                    r_scale = 1
            else:  # assume same scaling as polylines
                r_scale = 1
            if self.metadata.software == "smartroot":
                r_scale = scale
                print("DataModel.scale_rsml() radius length scale (smartroot)", r_scale)
            for i in range (0, len(self.radii)):
                self.radii[i] *= r_scale
            # creation times
            cts_scale = 1.  # assume it is in days
            if self.tagnames[1]:
                if self.tagnames[1] in self.metadata.properties:
                    cts_scale = self.scales_[self.metadata.properties[self.tagnames[1]].unit]
                    print("DataModel.scale_rsml() temporal scale", r_scale)
            for i in range (0, len(self.cts)):
                self.cts[i] *= cts_scale
            # polylines
            for i in range(0, len(self.polylines)):
                for j in range(0, len(self.polylines[i])):
                    for k in range(0, 3):
                        self.polylines[i][j][k] *= scale
        else:
            tkinter.messagebox.showwarning("Warning", "Open RSML file first")

    def convert_to_analyser(self):
        """ 
        converts the polylines to a SegmentAnalyser and a MappedSegments object, and stores max_ct 
        """
        if "parent-poly" in self.properties:
            nodes, segs = rsml_reader.get_segments(self.polylines, self.properties)  # fetch nodes and segments
            segRadii = np.zeros((segs.shape[0], 1))  # convert to paramter per segment
            segCTs = np.zeros((segs.shape[0], 1))
            subTypes = np.zeros((segs.shape[0], 1))
            if np.isnan(self.cts[0]):  # nan indicates creation times not given in rsml
                self.cts = np.zeros((len(self.cts),))
            for i, s in enumerate(segs):
                segRadii[i] = self.radii[s[1] - 1]  # seg to node index
                segCTs[i] = self.cts[s[1] - 1]
                subTypes[i] = self.types[s[1] - 1]
            if np.isnan(subTypes[0]):
                subTypes = np.ones((len(segs),), dtype=np.int64)
            self.max_ct = np.max(segCTs)
            segs_ = [pb.Vector2i(s[0], s[1]) for s in segs]  # convert to CPlantBox types
            nodes_ = [pb.Vector3d(n[0], n[1], n[2]) for n in nodes]
            self.analyser = pb.SegmentAnalyser(nodes_, segs_, segCTs, segRadii)
            minb = self.analyser.getMinBounds()
            maxb = self.analyser.getMaxBounds()
            if maxb.z > 0 and minb.z > 0:  # image coordinates in px often start in lower left corner
                print("DataModel.convert_to_analyser() assuming image coordinates, y-centered and z-flipped ")
                for i in range(len(nodes_)):
                    nodes_[i].z = -nodes_[i].z
                    nodes_[i].y -= maxb.y / 2  # center
                self.analyser = pb.SegmentAnalyser(nodes_, segs_, segCTs, segRadii)  # redefine with shifted nodes
            self.analyser.addData("subType", subTypes)
            self.mapped_segments = pb.MappedSegments(self.analyser.nodes, np.array(self.cts), segs_, np.array(segRadii), np.array(subTypes))
        else:
            print(properties.keys())  # 'parent-poly' should have been added by the reader
            tkinter.messagebox.showwarning("Warning", "'parent-poly' property is missing, cannot create SegmentAnalyser")

