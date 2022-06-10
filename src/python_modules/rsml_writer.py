import xml.etree.ElementTree as ET
from xml.dom import minidom

import datetime
from scipy import sparse
import numpy as np

""" 
RSML Writer, by Daniel Leitner (2019-2021) 
"""


class Property:
    """RSML property tag helper class, for property and function meta data description """

    def __init__(self, label:str, type_:str, unit:str, data:list):
        """
        Args:
        label(str): name of the property 
        type_(str): type ("int" or "float")
        unit(str): units (e.g. "cm", "days")
        data(list): property data per branch number
        """
        self.label = label
        self.type = type_
        self.unit = unit
        self.data = data

    def write_property(self, basetag:ET.Element):
        """ adds the rsml property tag to the basetag 
        
        Args:
        basetag(ET.Element): base tag in the xml element tree
        """
        base = ET.SubElement(basetag, "property-definition")
        ET.SubElement(base, "label").text = self.label
        ET.SubElement(base, "type").text = self.type
        ET.SubElement(base, "unit").text = self.unit


class Metadata:
    """RSML metadata tag helper class, including suggested default values, 
    all attributes are public, modify at will"""

    def __init__(self):
        self.unit = "cm"
        self.resolution = "1"
        self.last_modified = datetime.date.today().strftime("%d-%B-%Y")
        self.software = "Daniel's friendly RSML converter"
        self.image_label = ""  # if available put source image name here
        self.time_sequence_label = "artifical"
        self.time_sequence_index = "1/1"
        self.time_sequence_unified = "true"
        self.properties = {}
        self.func_names = []
        self.scale_to_cm = 1.
        self.set_scale_()

    def set_scale_(self):
        if self.unit == "cm":
            to_cm = 1.
        elif self.unit == "pixel" or self.unit == "px":
            to_cm = 1.
        elif self.unit == "inch":
            to_cm = 2.54
        elif self.unit == "m":
            to_cm = 100.
        elif self.unit == "mm":
            to_cm = 0.1
        else:
            print(self.unit)
            raise "Metadata.set_scale_: do not know unit " + self.unit
        self.scale_to_cm = to_cm / float(self.resolution)

    def add_property(self, prop:Property):
        """ add a property description 
        
        Args:
        prop(Property): the meta-data desription of a property
        """
        self.properties[prop.label] = prop

    def set_fun_names(self, names:list):
        """ adds function names to the meta data
        
        Args:
        names(list) : list of strings containing the function names
        """
        self.func_names = names

    def write_meta(self, basetag:ET.Element):
        """ adds the rsml metadata tag to the basetag
        
            Args:
            basetag(ET.Element): base tag in the xml element tree
        """
        metadata = ET.SubElement(basetag, "metadata")
        ET.SubElement(metadata, "version").text = "1"
        ET.SubElement(metadata, "unit").text = self.unit
        ET.SubElement(metadata, "resolution").text = self.resolution
        ET.SubElement(metadata, "last-modified").text = self.last_modified
        ET.SubElement(metadata, "software").text = self.software
        if self.image_label != "":
            image = ET.SubElement(metadata, "image")
            ET.SubElement(image, "label").text = self.image_label
        # Properties
        prop_defs = ET.SubElement(metadata, "property-definitions")
        for p in self.properties.values():
            p.write_property(prop_defs)
        time_seq = ET.SubElement(metadata, "time-sequence")
        # Time Sequence
        ET.SubElement(time_seq, "label").text = self.time_sequence_label
        ET.SubElement(time_seq, "index").text = self.time_sequence_index
        ET.SubElement(time_seq, "unified").text = self.time_sequence_unified

    def read_meta(self, metadata_tag: ET.Element):
        """ reads from RSML metadata tag         
            Currently, only 'unit' and 'resolution' are read (TODO)        
        """
        for unit in metadata_tag.iterfind('unit'):
            self.unit = unit.text
        self.resolution = 1  # default
        for res in metadata_tag.iterfind('resolution'):
            self.resolution = res.text
        self.set_scale_()
        self.software = "unknown"
        for s in metadata_tag.iterfind('software'):
            self.software = s.text
        for pd in metadata_tag.iterfind('property-definitions'):
            for p in pd:
                label = ""
                for l in p.iterfind('label'):
                    label = l.text
                type_ = ""
                for v in p.iterfind('type'):
                    type_ = v.text
                unit = "1"
                for u in p.iterfind('unit'):
                    unit = u.text
                if label:
                    self.add_property(Property(label, type_, unit, []))


class LinkedPolylines:
    """ LinkedPolylines, helper class for writing rsml files, used by follow_, which is used by segs2polylines
    
    Attributes:
        branchnumer(int): the branch number of this polyline
        polylines(list of int): holds the node indices of this polyline
        lateral(list of LinkedPolyLines) references to the laterals roots, last element is an empty list
        data(list of list of float): optionally, holds additional data per root for rsml properties       
        base_root (bool): if it is a base root, or not
        bc(static int): branch counter
        metadata(static Metadata): meta data of the root system, to look up rsml property and function names
    """

    bc = 0  # static member variable
    metadata = []  # static member variable

    def __init__(self):
        self.branchnumber = -1
        self.parent_node = -1  # parent node index!
        self.polyline = []
        self.laterals = []
        self.base_root = False

    @staticmethod
    def set_metadata(meta):
        """ Sets the meta data of the root system, to look up property and function names 
        
        Args:
        meta(Metadata): the rsml metadata 
        """
        LinkedPolylines.metadata = meta

    def write_root(self, basetag, nodes, nodedata, **kwargs):
        """ Adds the rsml root tag to the basetag
        """
        meta = LinkedPolylines.metadata  # rename
        LinkedPolylines.bc += 1
        if kwargs.get("Renumber", True):
            root = ET.SubElement(basetag, "root", ID = str(LinkedPolylines.bc))
        else:
            root = ET.SubElement(basetag, "root", ID = str(self.branchnumber))
        # Geometry
        geom = ET.SubElement(root, "geometry")
        pl = ET.SubElement(geom, "polyline")
        for p in self.polyline:
            x_ = nodes[p, 0]
            y_ = nodes[p, 1]
            z_ = nodes[p, 2]
            pl.append(ET.Element("Point", dict(x = str(x_), y = str(y_), z = str(z_))))
        # Properties
        properties = ET.SubElement(root, "properties")  # defined by root
        for p in meta.properties.values():
            if p.data:
                properties.append(ET.Element(p.label, dict(value = str(p.data[self.branchnumber - 1]))))
        if self.base_root:
            properties.append(ET.Element("parent-node", dict(value = str(self.parent_node))))
        else:
            properties.append(ET.Element("parent-node", dict(value = str(self.parent_node - 1))))
        for r in self.laterals:
            if r:  # ends with an empty list
                r.write_root(root, nodes, nodedata, Renumber = kwargs.get("Renumber", True))
        # Functions
        funcs = ET.SubElement(root, "functions")
        for i, n in enumerate(meta.func_names):
            fun = ET.SubElement(funcs, "functions", dict(domain = "polyline", name = n))  # defined by node
            for p in self.polyline:
                fun.append(ET.Element("sample", dict(value = str(nodedata[i][p]))))


def follow_(i0:int, i:int, i_old, A) -> LinkedPolylines:
    """ Recursively follows the roots
    
    Args:
        i0 (int): parent node index
        i (int): initial node  
        i_old (int) : parent node
        A (sparse matrix): adjacency matrix with branch ids, values are branch number (>0!) or root order (>0!) etc.
            
    Returns: 
        LinkedPolyline: representing the root system
    """
    _, j_ = A[i,:].nonzero()

    newlines = LinkedPolylines()
    newlines.polyline = [i]
    newlines.parent_node = i0

    if i_old > -1:  # in case it is not the first node
        newlines.branchnumber = A[i_old, i]
    else:
        newlines.branchnumber = A[i, j_[0]]  # pick first segment

    if len(j_) > 0:  # if there is at least another node

        first = True
        done = False
        while not done:
            done = True
            for j in j_:  # we process all neighbouring nodes
                if  A[i, j] == newlines.branchnumber and first:  # follow this branch number (maximal one node can have the same branchnumber)
                    newlines.polyline.append(j)
                    done = False  # if no other node has the same branch number, we quit
                    first = False
                else:  # laterals are created into all other directions
                    if  A[i, j] == newlines.branchnumber:
                        print("rsml_writer.follow_() warning: path is not unique, e.g. if root and lateral have the same branch id (order, branchnumber, ...)")
                    newlines.laterals.append(follow_(len(newlines.polyline) - 1, j, i, A))  # follow lateral
            if not done:
                i = newlines.polyline[-1]  # jump to next node
                _, j_ = A[i,:].nonzero()
                first = True

    return newlines


def segs2polylines(axes:list, segs:list, segdata:list) -> list:
    """ Converts a root system represented by nodes and segments into a linked polyline representation.
    
    Args:
        axes (ĺist of int): indices of the starting nodes
        segs (list of list of int) segments represented by two node indices 
        segdata (list of int) branch number or root order per segment (both >0!)
        
    Returns: 
        A list of LinkedPolyline: representing the root system axes    
    """
    n = max(max(segs[:, 0]), max(segs[:, 1])) + 1
    A = sparse.csr_matrix((segdata, (segs[:, 0], segs[:, 1])), dtype = np.int64, shape = (n, n))
    polylines = []
    for a in axes:
        polylines.append(follow_(-1, a, -1, A))
        polylines[-1].base_root = True  # follow_ adds subtree, main axis is last (?)
    return polylines


def write_rsml(name:str, axes:list, segs:list, segdata:list, nodes:list, nodedata:list, meta:Metadata, **kwargs):
    """ RSML writer
    
    Args:
        name(str): file name
        axes(list of int): indices of the starting nodes
        segs(list of list of int) segments: represented by two node indices 
        segdata(list of int): branch id, e.g. branch number or root order per segment (both >0!), needed to reconstruct the polylines from the segments
        nodes(list of list of float) nodes of the root system
        nodedata (list of list of float) data attached to the nodes, stored as functions
        meta(Metadata): Additional information, e.g. on rsml properties and functions
        
    kwargs: "Renumber" = True recalculates branch numbers
    """
    rsml = ET.Element("rsml")
    # Metadata
    meta.write_meta(rsml)
    # Scene
    LinkedPolylines.set_metadata(meta)
    polylines = segs2polylines(axes, segs, segdata)

    scene = ET.SubElement(rsml, "scene")
    plant = ET.SubElement(scene, "plant")
    for root in polylines:
        if root:
            root.write_root(plant, nodes, nodedata, Renumber = kwargs.get("Renumber", True))
    # Write (pretty)
    xmlstr = minidom.parseString(ET.tostring(rsml)).toprettyxml(indent = "   ", encoding = 'UTF-8')
    with open(name, "w") as f:
        f.write(xmlstr.decode('UTF-8'))
        f.close()


if __name__ == "__main__":

    axis = [0]  # initial node
    nodes = np.array([ [0.00, 0.00, -3.00], [-0.00, -0.01, -3.48], [-0.85, 0.48, -3.71], [-1.69, 0.99, -3.90], [-2.58, 1.32, -4.21], [-3.48, 1.67, -4.49], [-4.38, 2.00, -4.77], [-5.24, 2.40, -5.09], [-6.08, 2.82, -5.42], [-6.93, 3.27, -5.69], [-6.96, 3.29, -5.70], [-0.00, 0.01, -3.97], [0.20, -0.95, -4.20], [0.43, -1.88, -4.49], [0.65, -2.81, -4.77], [0.84, -3.75, -5.06], [1.04, -4.70, -5.31], [1.27, -5.64, -5.54], [1.43, -6.58, -5.84], [1.48, -6.91, -5.94], [-0.01, 0.03, -4.45], [0.75, 0.68, -4.48], [1.52, 1.32, -4.50], [2.30, 1.94, -4.46], [3.07, 2.58, -4.41], [3.88, 3.16, -4.46], [4.73, 3.69, -4.50], [5.34, 4.05, -4.53], [-0.03, 0.06, -4.97], [-0.73, 0.63, -5.40], [-1.46, 1.20, -5.79], [-2.15, 1.80, -6.18], [-2.76, 2.48, -6.59], [-3.16, 3.16, -7.21], [-3.63, 3.90, -7.64], [-0.06, 0.07, -5.42], [0.07, 1.04, -5.61], [0.23, 2.00, -5.82], [0.46, 2.93, -6.11], [0.72, 3.85, -6.40], [0.99, 4.77, -6.69], [1.06, 5.11, -6.78], [-0.10, 0.08, -5.96], [0.35, 0.90, -6.31], [0.80, 1.72, -6.65], [1.23, 2.55, -7.01], [1.62, 3.41, -7.34], [1.70, 3.58, -7.41], [-0.14, 0.10, -6.46], [0.05, -0.81, -6.83], [0.29, -1.72, -7.18], [0.43, -2.62, -7.59], [0.50, -3.12, -7.84], [-0.19, 0.12, -7.02], [0.32, 0.98, -7.03], [0.82, 1.84, -6.95], [1.07, 2.23, -6.88], [-0.23, 0.17, -7.56], [0.68, -0.18, -7.80], [1.12, -0.35, -7.93], [-0.24, 0.21, -8.14], [-0.01, 0.36, -8.22], [-0.25, 0.24, -8.69], [-0.24, 0.25, -9.25], [-0.24, 0.26, -9.71], [-0.26, 0.26, -10.09], [-0.28, 0.25, -10.57], [-0.26, 0.24, -11.05], [-0.26, 0.21, -11.58], [-0.25, 0.19, -12.06], [-0.25, 0.17, -12.55], [-0.23, 0.15, -13.00], [-0.23, 0.12, -13.46], [-0.22, 0.12, -13.99], [-0.21, 0.15, -14.54], [-0.19, 0.20, -15.07], [-0.17, 0.24, -15.60], [-0.17, 0.31, -16.17], [-0.15, 0.36, -16.64], [-0.12, 0.43, -17.18], [-0.11, 0.48, -17.70], [-0.07, 0.52, -18.23], [-0.06, 0.53, -18.53] ])
    segs = np.array([ [0, 1], [1, 11], [11, 20], [20, 28], [28, 35], [35, 42], [42, 48], [48, 53], [53, 57], [57, 60], [60, 62], [62, 63], [63, 64], [64, 65], [65, 66], [66, 67], [67, 68], [68, 69], [69, 70], [70, 71], [71, 72], [72, 73], [73, 74], [74, 75], [75, 76], [76, 77], [77, 78], [78, 79], [79, 80], [80, 81], [81, 82], [1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7], [7, 8], [8, 9], [9, 10], [11, 12], [12, 13], [13, 14], [14, 15], [15, 16], [16, 17], [17, 18], [18, 19], [20, 21], [21, 22], [22, 23], [23, 24], [24, 25], [25, 26], [26, 27], [28, 29], [29, 30], [30, 31], [31, 32], [32, 33], [33, 34], [35, 36], [36, 37], [37, 38], [38, 39], [39, 40], [40, 41], [42, 43], [43, 44], [44, 45], [45, 46], [46, 47], [48, 49], [49, 50], [50, 51], [51, 52], [53, 54], [54, 55], [55, 56], [57, 58], [58, 59], [60, 61] ])
    age = np.array([ 8.00, 7.76, 7.52, 7.29, 7.03, 6.80, 6.53, 6.28, 6.00, 5.73, 5.43, 5.16, 4.87, 4.64, 4.44, 4.20, 3.95, 3.68, 3.43, 3.17, 2.94, 2.70, 2.42, 2.14, 1.85, 1.57, 1.27, 1.02, 0.73, 0.45, 0.16, -0.00, 2.29, 2.02, 1.74, 1.43, 1.11, 0.77, 0.41, 0.01, 0.00, 2.04, 1.77, 1.48, 1.18, 0.85, 0.51, 0.14, 0.00, 1.78, 1.51, 1.22, 0.92, 0.60, 0.26, 0.00, 1.51, 1.24, 0.96, 0.66, 0.34, 0.00, 1.28, 1.01, 0.73, 0.43, 0.12, 0.00, 0.97, 0.69, 0.39, 0.07, 0.00, 0.72, 0.45, 0.17, 0.00, 0.41, 0.14, 0.00, 0.13, 0.00, 0.00 ])
    types = np.array([ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2 ])
    a_tap = 0.2  # tap root radius (cm)
    a_lateral = 0.1  # lateral root radius (cm)

    meta = Metadata()
    meta.add_property(Property("radius", "float", "cm", [a_tap, a_lateral]))
    meta.set_fun_names(["age"])
    write_rsml("test.xml", axis, segs, types, nodes, [age], meta)
    print("done")
