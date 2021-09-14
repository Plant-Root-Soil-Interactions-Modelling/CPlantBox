import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt
import numpy as np
import math

""" 
    RSML reader by Daniel Leitner 2019
    usage:
    polylines, properties, functions = read_rsml("RootSystem.rsml")
"""


def parse_rsml(organ :ET, polylines :list, properties :dict, functions :dict, parent :int) -> (list, dict, dict):
    """Recursivly parses the rsml file 
    """
    for poly in organ.iterfind('geometry'):  # only one
        polyline = []
        for p in poly[0]:  # 0 is the polyline
            n = p.attrib
            newnode = [ float(n['x']), float(n['y']), float(n['z']) ]
            polyline.append(newnode)
        polylines.append(polyline)
        properties.setdefault("parent-poly", []).append(parent)

    for prop in organ.iterfind('properties'):
        for p in prop:  # i.e legnth, type, etc..
            properties.setdefault(str(p.tag), []).append(float(p.attrib['value']))

    for funcs in organ.iterfind('functions'):
        for fun in funcs:
            samples = [ ]
            for sample in fun.iterfind('sample'):
                samples.append(float(sample.attrib['value']))
            functions.setdefault(str(fun.attrib['name']), []).append(samples)

    pi = len(polylines) - 1
    for elem in organ.iterfind('root'):  # and all laterals
        polylines, properties, functions = parse_rsml(elem, polylines, properties, functions, pi)

    return polylines, properties, functions


def read_rsml(name :str) -> (list, dict, dict):
    """Parses the RSML file into: 
    (1) a (flat) list of polylines, with one polyline per root
    (2) a dictionary of properties one per root, 
    adds "parent_poly" holding the index of the parent root in the list of polylines
    (3) a dictionary of functions     
    """
    root = ET.parse(name).getroot()
    plant = root[1][0]
    polylines = []
    properties = { }
    functions = { }
    for elem in plant.iterfind('root'):
        (polylines, properties, functions) = parse_rsml(elem, polylines, properties, functions, -1)

    return polylines, properties, functions


def get_segments(polylines :list, props :dict) -> (list, list):
    """ Converts the polylines to a list of nodes 
    and an index list of line segments 
    """
    nodes, offset, segs = [], [], []
    offset.append(0)  # global node index at each polyline
    for p in polylines:
        for n in p:
            nodes.append(n)
        offset.append(offset[-1] + len(p))
    for i, p in enumerate(polylines):
        ni = props["parent-node"][i]
        pi = props["parent-poly"][i]
        if (pi >= 0):
            segs.append([offset[pi] + ni, offset[i]])
        for j in range(0, len(p) - 1):
            segs.append([offset[i] + j, offset[i] + j + 1])
    return nodes, segs


def get_parameter(polylines :list, funcs :dict, props :dict) -> (list, list, list):
    """ Copies radii and creation times, one value per segment 
    """
    fdiam = funcs["diameter"]
    fet = funcs["emergence_time"]
    ptype = props["type"]
    radii, cts, types = [], [], []
    for i, p in enumerate(polylines):
        for j in range(0, len(p)):
            radii.append(fdiam[i][j] / 2)
            cts.append(fet[i][j])
            types.append(ptype[i])
    return radii, cts, types


def plot_rsml(polylines :list, prop : list):
    """Plots the polylines in y-z axis
    """
    cmap = plt.get_cmap("hsv", 256)
    newcolors = cmap(np.linspace(0, 1, math.ceil(max(prop)) + 1))
    for i, pl in enumerate(polylines):
        nodes = np.array(pl)
        plt.plot(nodes[:, 1], nodes[:, 2], color = newcolors[int(prop[i]), :])  # y,z plot / (len(polylines) - 1)
    plt.axis('equal')
    plt.show()


def plot_segs(nodes, segs):
    """Plots the segments in y-z axis
    """
    for s in segs:
        plt.plot([nodes[s[0], 1], nodes[s[1], 1]], [nodes[s[0], 2], nodes[s[1], 2]], "r")
    plt.axis('equal')
    plt.show()


if __name__ == '__main__':
    polylines, properties, functions = read_rsml("organism.rsml")
    print("Properties:")
    for key, v in properties.items() :
        print("\t", key, len(properties[key]))
    print("Functions:")
    for key, v in functions.items() :
        print("\t", key, len(functions[key]))
    nodes, segs = get_segments(polylines, properties)
    radii, cts, types = get_parameter(polylines, functions, properties)

    # plot_rsml(polylines, properties["parent-node"])  # properties["parent-node"]
    nodes = np.array(nodes)
    segs = np.array(segs, dtype = np.int64)
    plot_segs(nodes, segs)

