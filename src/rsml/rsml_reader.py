import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from scipy.spatial import distance

from rsml.rsml_writer import Metadata

""" 
RSML Reader, by Daniel Leitner (2019-2021), see also rsml_data.py
"""


def parse_rsml_(organ:ET, polylines:list, properties:dict, functions:dict, parent:int) -> (list, dict, dict):
    """ Recursivly parses the rsml file, used by read_rsml """
    for poly in organ.iterfind('geometry'):  # only one
        polyline = []
        for p in poly[0]:  # 0 is the polyline
            n = p.attrib
            if 'z' in n:
                newnode = [ float(n['x']), float(n['y']), float(n['z']) ]
            else:
                newnode = [ 0., float(n['x']), float(n['y']) ]
            polyline.append(newnode)
        polylines.append(polyline)
        properties.setdefault("parent-poly", []).append(parent)

    for prop in organ.iterfind('properties'):
        for p in prop:  # i.e legnth, type, etc..
            if 'value' in p.attrib:
                try:
                    properties.setdefault(str(p.tag), []).append(float(p.attrib['value']))
                except:
                    properties.setdefault(str(p.tag), []).append(np.nan)
            else:
                properties.setdefault(str(p.tag), []).append(float(p.text))
    for funcs in organ.iterfind('functions'):
        for fun in funcs:
            samples = [ ]
            for sample in fun.iterfind('sample'):
                if 'value' in sample.attrib:
                    samples.append(float(sample.attrib['value']))
                else:
                    samples.append(float(sample.text))
            functions.setdefault(str(fun.attrib['name']), []).append(samples)

    pi = len(polylines) - 1
    for elem in organ.iterfind('root'):  # and all laterals
        polylines, properties, functions = parse_rsml_(elem, polylines, properties, functions, pi)

    return polylines, properties, functions


def read_rsml(name:str) -> (list, dict, dict, Metadata):
    """ Parses the RSML file into:

    Args:
    name(str): file name of the rsml file

    Returns: 
    (list, dict, dict, Metadata):
    (1) a (flat) list of polylines, with one polyline per root
    (2) a dictionary of properties, one per root, adds "parent_poly" holding the index of the parent root in the list of polylines
    (3) a dictionary of functions  
    (4) metadata   
    """
    root = ET.parse(name).getroot()

    metadata = Metadata()
    metadata.read_meta(root[0])

    plant = root[1][0]
    polylines = []
    properties = { }
    functions = { }
    for elem in plant.iterfind('root'):
        (polylines, properties, functions) = parse_rsml_(elem, polylines, properties, functions, -1)

    if metadata.software == "RSWMS":
         del properties['parent-node']  # they are something unexpected

    return polylines, properties, functions, metadata


def artificial_shoot(polylines, properties, functions):
    """ inserts an artificial shoot, with functions and properties of the the first polyline     
    """
    polylines.insert(0, [[0, 0, -0.1], [0, 0, -2.9]])  # artificial shoot
    for key, v in properties.items():
        properties[key].insert(0, properties[key][0])
    for key, v in functions.items():
        functions[key].insert(0, [functions[key][0][0], functions[key][0][1]])
    for i, p in enumerate(polylines):  # add one to the parent poly indices
        properties["parent-poly"][i] += 1
    for i, p in enumerate(polylines):
        if properties["parent-poly"][i] == 0:
            properties["parent-node"][i] = 1
    properties["parent-poly"][0] = -1
    properties["parent-node"][0] = -1


def get_segments(polylines:list, props:dict) -> (list, list):
    """ Converts the polylines to a list of nodes and an index list of line segments
        
    Args:
    polylines(list): flat list of polylines, one polyline per root
    props(dict): dictionary of properties, one per root, must contain "parent-node", (and "parent-poly" that was added by read_rsml)
    
    Returns: 
    (list, list): 
    (1) list of nodes
    (2) list of two integer node indices for each line segment 
    """
    nodes, offset, segs = [], [], []
    offset.append(0)  # global node index at each polyline
    if not "parent-node" in props:
        print("rsml_reader.get_segments(): no 'parent-node' tag found using nearest node")
        add_parent_nodes(polylines, props)
    for p in polylines:
        for n in p:
            nodes.append(n)
        offset.append(offset[-1] + len(p))
    for i, p in enumerate(polylines):
        pi = props["parent-poly"][i]
        if (pi >= 0):  # not a base root
            ni = int(props["parent-node"][i])
            # print(i, pi, len(polylines[pi]), ni)
            # print(polylines[i])
            # print(polylines[pi])
            # print(props["subType"][i])
            # print(props["subType"][pi])
            assert ni < len(polylines[pi]), "parent node index exceeds number of parent nodes"
            segs.append([offset[pi] + ni, offset[i]])
        for j in range(0, len(p) - 1):
            segs.append([offset[i] + j, offset[i] + j + 1])
    return np.array(nodes), np.array(segs, dtype = np.int64)


def add_parent_nodes(polylines, props):
    """ adds the parent-node property, by minimizing the distance to the first point of the lateral """
    props["parent-node"] = [None] * len(polylines)  # empty list
    for i, p in enumerate(polylines):
        x = np.array([p[0]])
        pi = props["parent-poly"][i]
        if (pi >= 0):
            # print(pi, len(polylines))
            y = np.array(polylines[pi])
            props["parent-node"][i] = np.argmin(distance.cdist(x, y))
        else:
            props["parent-node"][i] = -1


def age_to_creationtime(age):
    """ ages per polyline as list of list """
    maxage = 0.
    for pl in age:
        for a in pl:
            maxage = max(maxage, a)
    return [[maxage - a for a in pl] for pl in age]


def get_root_parameters(polylines:list, funcs:dict, props:dict) -> (list, list, list):
    """ Copies radius, creation times, and types one value per root, if type or order is not found, root order is calculated 
        
        Args:
        polylines(list): flat list of polylines, one polyline per root
        funcs(dict): dictionary of functions
        props(dict): dictionary of properties     
        
        Returns:
        radius, creation time, type or order (per root), selected tag names, radii is single valued, type is single valued
        
        different rsml writer call parameters different names, get_parameter expect the ones in the list below
        functions are checked first, then properties, if not found NaN values are set.       
        
    """
    radius_names = ["radius", "radii", "Radius", "Radii", "radius [m]"]  # add more, where needed
    diam_names = ["diameter", "diameters", "diam", "Diameter", "Diameters", "rootDiameter"]
    type_names = ["type", "types", "subType", "subTypes", "order", "orders"]
    ct_names = ["creation_time", "creationTime", "emergence_time", "emergenceTime", "node_creation_time", "nodeCreationTime"]
    age_names = ["age", "Age", "age [d]"]

    tag_names = []
    radii = None
    for n in radius_names:  # search for radius
        if n in funcs:  # in functions
            tag_names.append(n)
            radius = True
            radii = funcs[n]
            radii_p = False
            break
        elif n in props:  # and in properties
            tag_names.append(n)
            radius = True
            radii = props[n]
            radii_p = True
            break
    if radii == None:  # nothing found yet
        for n in diam_names:  # search for radiieter
            if n in funcs:  # in functions
                tag_names.append(n)
                radius = False
                radii = funcs[n]
                radii_p = False
                break
            elif n in props:  # and in properties
                tag_names.append(n)
                radius = False
                radii = props[n]
                radii_p = True
                break
    if radii == None:  # nothing found
        tag_names.append("")

    et = None
    for n in ct_names:  # search for creation times
        if n in funcs:  # only in functions
            tag_names.append(n)
            et = funcs[n]
            break
    if et == None:  # nothing found yet
        for n in age_names:  # search for creation times
            if n in funcs:  # only in functions
                tag_names.append(n)
                et = age_to_creationtime(funcs[n])
                break
    if et == None:  # nothing found
        tag_names.append("")

    type_ = None
    for n in type_names:  # search for types
        if n in funcs:  # in functions
            tag_names.append(n)
            type_ = funcs[n]
            type_p = False
            break
        elif n in props:  # and in properties
            tag_names.append(n)
            type_ = props[n]
            type_p = True
            break
    if type_ == None:  # nothing found
        type_ = get_root_orders(props)
        type_p = True
        props["order"] = type_
        tag_names.append("order")

    if radius == False:
        for i, p in enumerate(polylines):
            for j in range(0, len(p)):
                if radii_p:
                    radii[i] = radii[i] / 2.
                else:
                    radii[i][j] = radii[i][j] / 2.

    return radii, et, type_, tag_names, radii_p, type_p


def get_parameter(polylines:list, funcs:dict, props:dict) -> (list, list, list):
    """ Copies radius, creation times, and types one value per node, if type or order is not found, root order is calculated.
        see also get_root_parameters
        
        Args:
        polylines(list): flat list of polylines, one polyline per root
        funcs(dict): dictionary of functions
        props(dict): dictionary of properties     
        
        Returns:
        radius, creation time, type or order (per node)                                    
    """

    radii_, et, type_, tag_names, radii_p, type_p = get_root_parameters(polylines, funcs, props)

    radii, types, cts = [], [], []
    for i, p in enumerate(polylines):
        for j in range(0, len(p)):
            if radii_ is not None:
                if radii_p:
                    radii.append(radii_[i])
                else:
                    radii.append(radii_[i][j])
            else:
                radii.append(np.nan)
            if type_ is not None:
                if type_p:
                    types.append(type_[i])
                else:
                    types.append(type_[i][j])
            else:
                types.append(np.nan)
            if et is not None:
                cts.append(et[i][j])
            else:
                cts.append(np.nan)

    return radii, cts, types, tag_names


def get_property(name:str, polylines:list, props:dict) -> (list):
    """ Copies a rsml property to one value per node 
        
        Args:
        name(str): name of the rsml properyt
        polylines(list): flat list of polylines, one polyline per root
        props(dict): dictionary of properties     
        
        Returns:
        the property converted to per node        
    """
    p = []
    for i, p in enumerate(polylines):
        for j in range(0, len(p)):
            p.append(props[name])
    return p


def order_(pp, i, c = 0):
    """ 
    recursively finds the root order, seee get_root_orders()
    """
    if pp[i] == -1:
        return c
    else:
        return order_(pp, pp[i], c + 1)


def get_root_orders(props):
    """ 
    calculates all root orders based on the 'parent-poly' tag 
    """
    pp = props["parent-poly"]
    orders = []
    for i in range(0, len(pp)):
        orders.append(order_(pp, i))
    return orders


def plot_rsml(polylines:list, prop:list):
    """Plots the polylines in y-z axis with colors given by a root property

    Args:
    polylines(list): flat list of polylines, one polyline per root 
    prop(list): a single property, list of scalar value, on per root 
    """
    f = matplotlib.colors.Normalize(vmin = min(prop), vmax = max(prop))
    cmap = plt.get_cmap("jet", 256)
    for i, pl in enumerate(polylines):
        nodes = np.array(pl)
        plt.plot(nodes[:, 1], nodes[:, 2], color = cmap(f(prop[i])))
    plt.axis('equal')
    plt.show()


def plot_multiple_rsml(polylines:list, prop:list, times):
    """Plots the polylines in y-z axis with colors given by a root property
    Args:
    polylines(list): flat list of polylines, one polyline per root 
    prop(list): a single property, list of scalar value, on per root 
    """
    n = len(polylines)
    fig, axes = plt.subplots(1, n, figsize = (15, 7))
    f = matplotlib.colors.Normalize(vmin = min(prop[-1]), vmax = max(prop[-1]))
    cmap = plt.get_cmap("jet", 256)
    for j in range(0, n):
        for i, pl in enumerate(polylines[j]):
            nodes = np.array(pl)
            if i == 0:
                axes[j].plot(nodes[:, 1], nodes[:, 2], "r")  # y,z plot
            else:
                axes[j].plot(nodes[:, 1], nodes[:, 2], color = cmap(f(prop[j][i])))  # y,z plot
        axes[j].set_xlim([-5, 5.])
        axes[j].set_ylim([-15., 0.])
        axes[j].set_title("{} days".format(times[j]))
    fig.tight_layout()


def plot_segs(nodes:list, segs:list, fun:list):
    """Plots the segments in y-z axis (rather slow)
    
    Args:
    nodes(list): list of nodes
    segs(list): list of two integer node indices for each line segment 
    fun(list): a single function, list of scalar value, on per segment, see TODO 
    """
    f = matplotlib.colors.Normalize(vmin = min(fun), vmax = max(fun))
    cmap = plt.get_cmap("jet", 256)
    print("Segments")
    for i, s in enumerate(segs):
        plt.plot([nodes[s[0], 1], nodes[s[1], 1]], [nodes[s[0], 2], nodes[s[1], 2]], color = cmap(f(fun[i])))
    plt.axis('equal')
    plt.show()


if __name__ == '__main__':

    # fname = "../../../dumux-rosi/grids/RootSystem.rsml"
    fname = "../../tutorial/examples/python/results/example_3c.rsml"  # run example3c_write.py first (in tutorial/examples/python/)

    polylines, properties, functions, _ = read_rsml(fname)
    artificial_shoot(polylines, properties, functions)  # for multiple base roots, add artificial root

    print("Properties:")
    for key, v in properties.items():
        print("\t", key, len(properties[key]))
    print("Functions:")
    for key, v in functions.items():
        print("\t", key, len(functions[key]))
    plot_rsml(polylines, properties["parent-node"])

    nodes, segs = get_segments(polylines, properties)
    nodes = np.array(nodes)
    segs = np.array(segs, dtype = np.int64)

    radii, cts, types, tag_names = get_parameter(polylines, functions, properties)
    plot_segs(nodes, segs, cts)  # slow
