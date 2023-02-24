import sys; sys.path.append("../../src/python_modules/"); sys.path.append("../../")

import plantbox as pb
import functional.xylem_flux as xylem_flux
import rsml.rsml_reader as rsml_reader
from rsml.rsml_data import RsmlData

import numpy as np


class ViewerDataModel(RsmlData):
    """ 
    Model in the sense of model view controller (MVC), stores most that is presented in the view
    
    * manages the reading of rsml DataModel.open_rsml (see RsmlData defined in src/python_modules/rsml_data.py)
    * can add an artifical shoot, see DataModel.add_artificial_shoot
    * can interpolate creation time  DataModel.add_creation_times
    """

    def __init__(self):
        super().__init__()
        self.analyser = None  # created by convert_to_xylem_flux_
        self.xylem_flux = None  # created by convert_to_xylem_flux_
        self.base_nodes = [0]  # base nodes indices (of roots or multiple plants)
        self.base_segs = [0]  # emerging segment indices from base nodes

    def open_rsml(self, fname, z_shift = False):
        """ see RsmlData.open_rsml() in src/python_modules/rsml_data.py                
        Additionally, creates an analyser (pb.SegmentAnalyser) and a xylem_flux (pb.XylemFluxPython) object 
        """
        RsmlData.open_rsml(self, fname, z_shift)
        if not self.tagnames[0]:
            print("ViewerDataModel.open_rsml: no radius tag found, set to 0.1 cm")
            c = 0  # node counter
            for i, pl in enumerate(self.polylines):
                c += 1
                for p in pl:
                    c += 1
            self.radii = np.ones((c,)) * 0.1  # cm
        self.convert_to_xylem_flux_()

    def convert_to_xylem_flux_(self):
        """  converts the polylines to a SegmentAnalyser and a MappedSegments object, and stores max_ct   
        
        uses:
        properties["parent-poly"], properties["parent-nodes"]
        radii, cts, types                    
        """
        nodes, segs = rsml_reader.get_segments(self.polylines, self.properties)  # fetch nodes and segments
        segRadii = np.zeros((segs.shape[0], 1))  # convert to paramter per segment
        segCTs = np.zeros((segs.shape[0], 1))
        subTypes = np.zeros((segs.shape[0], 1))
        for i, s in enumerate(segs):
            segRadii[i] = self.radii[s[1]]  # seg to node index
            segCTs[i] = self.cts[s[1]]
            subTypes[i] = self.types[s[1]]
#             st = subTypes[i]  # had coded values...
#             if st == 0:
#                 segRadii[i] = 0.055  # seg to node index
#             elif st == 1:
#                 segRadii[i] = 0.03
#             elif st == 2:
#                 segRadii[i] = 0.02
        if np.isnan(subTypes[0]):
            subTypes = np.ones((len(segs),), dtype = np.int64)
        segs_ = [pb.Vector2i(s[0], s[1]) for s in segs]  # convert to CPlantBox types
        nodes_ = [pb.Vector3d(n[0], n[1], n[2]) for n in nodes]
        self.analyser = pb.SegmentAnalyser(nodes_, segs_, segCTs, segRadii)
        self.analyser.addData("subType", subTypes)
        ms = pb.MappedSegments(self.analyser.nodes, np.array(self.cts), segs_, np.array(segRadii), np.array(subTypes))
        self.xylem_flux = xylem_flux.XylemFluxPython(ms)
        self.base_nodes = self.get_base_node_indices_()
        self.xylem_flux.neumann_ind = self.base_nodes  # needed for suf
        self.xylem_flux.dirichlet_ind = self.base_nodes  # needed for krs
        self.base_segs = self.xylem_flux.find_base_segments()

    def get_base_node_indices_(self):
        """ get all node indices of base roots """
        c = 0
        bni = []
        for i, p in enumerate(self.polylines):
            if self.properties['parent-poly'][i] == -1:
                bni.append(c)
            c += len(p)  # next first node index
        return bni

    def add_artificial_shoot(self):
        """ adds a 1 cm shoot element, connecting all base roots
        the type for looking up conductivities is set to 10 
        """
        nodes = self.analyser.nodes
        bni = self.base_nodes
        print("DataModel.add_artificial_shoot() base node indices are", bni)
        mid = np.zeros((3,))
        for i in bni:
            mid += np.array([nodes[i].x, nodes[i].y, nodes[i].z])
        mid /= len(bni)
        print("DataModel.add_artificial_shoot() mid point is", mid)

        rsml_reader.artificial_shoot(self.polylines, self.properties, self.functions)  # append artifial shoot (default values)
        radii, cts, types, tagnames = rsml_reader.get_parameter(self.polylines, self.functions, self.properties)  # paramter per node
        # change default values from artificial shoot
        collar = mid.copy()
        mid[2] += 0.1  # shift a bit up
        collar[2] += 1.1  # 1 cm shoot length
        self.polylines[0][0] = collar
        self.polylines[0][1] = mid
        self.set_selected(radii, cts, types, tagnames)
        self.scale_selected_()
        self.radii[0] = 0.1  # cm
        self.radii[1] = 0.1  # cm
        self.types[0] = 10
        self.types[1] = 10
        # print("after ADD SHOOT")
        # print(self.polylines)
        # print(self.properties["parent-poly"])
        # print(self.properties["parent-node"])
        # print("radii", radii)
        # print("cts", cts)
        # print("types", types)
        self.convert_to_xylem_flux_()
#         print("mid ", str(self.analyser.nodes[1]), " cm")
#         print("collar ", str(self.analyser.nodes[0]), " cm")
#         print("seg 0", str(self.analyser.segments[0]))
#         print("radius", str(self.analyser.data["radius"][0]))
#         print("type", str(self.analyser.data["subType"][0]))
#         print("seg 1", str(self.analyser.segments[1]))
#         print("radius", str(self.analyser.data["radius"][1]))
#         print("type", str(self.analyser.data["subType"][1]))

    def add_creation_times(self):
        """ lineary interpolates creation times assuming a lateral delay time of one day """
        pl_ = self.polylines  # rename
        self.functions["creation_time"] = [None] * len(pl_)
        for i, pl in enumerate(pl_):
            if self.properties["parent-poly"][i] == -1:
                self.functions["creation_time"][i] = np.zeros((len(pl,)))
                if not (i == 0 and len(pl) == 2):  # not artifical shoot, else [0,0] is fine
                    ct = self.functions["creation_time"][i]  # rename
                    lt = self.get_length_(pl)
                    l = 0
                    ct[0] = 0
                    for j in range(0, len(pl) - 1):
                        l += np.linalg.norm(np.array(pl[j + 1]) - np.array(pl[j]))
                        ct[j + 1] = self.max_ct * (l / lt)
        for i, pl in enumerate(pl_):
            if self.functions["creation_time"][i] is None:
                self.add_creation_times_(i)
        radii, cts, types, tagnames = rsml_reader.get_parameter(self.polylines, self.functions, self.properties)  # paramter per node
        self.set_selected(radii, cts, types, tagnames)
        self.scale_selected_()
        self.convert_to_xylem_flux_()

    def add_creation_times_(self, i):
        """ recursive funtion called by add_creation_times, 
        interpolates a single polyline assuming a lateral delay time of one day
        """
        ppi = self.properties["parent-poly"][i]
        if self.functions["creation_time"][ppi] is None:
            print(ppi, self.functions["creation_time"][ppi])
            self.add_creation_times_(ppi)  # label parent first
        pl = self.polylines[i]
        self.functions["creation_time"][i] = np.zeros((len(pl),))
        ct = self.functions["creation_time"][i]  # rename
        lt = self.get_length_(pl)
        pni = self.properties["parent-node"][i]
        # print("Root", i, "Parent", ppi, "node index", pni)
        ct[0] = self.functions["creation_time"][ppi][pni]  #  + 1  # 1 day
        l = 0.
        for j in range(0, len(pl) - 1):
            l += np.linalg.norm(np.array(pl[j + 1]) - np.array(pl[j]))
            ct[j + 1] = ct[0] + (self.max_ct - ct[0]) * (l / lt)

    def get_length_(self, pl):
        """ polyline length """
        lt = 0
        for j in range(0, len(pl) - 1):
            lt += np.linalg.norm(np.array(pl[j + 1]) - np.array(pl[j]))
        return lt
