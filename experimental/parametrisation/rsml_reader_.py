import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt
import matplotlib
import numpy as np

""" 
RSML Reader, by Daniel Leitner (2019) 
"""


def parse_rsml_(organ:ET, polylines:list, properties:dict, functions:dict, parent:int) -> (list, dict, dict):
    """ Recursivly parses the rsml file, used by read_rsml """
    for poly in organ.iterfind('geometry'):  # only one
        polyline = []
        for p in poly[0]:  # 0 is the polyline
            n = p.attrib
            if len(n) == 3:  # RSML is 3D
                newnode = [ float(n['x']), float(n['y']), float(n['z']) ]
            else:  # RSML is 2D
                newnode = [ float(n['x']) / 116.93, 0.0, -1 * float(n['y']) / 116.93 + 10 ]  # 116.93 is resolution in pixel/cm TODO: get this info out of RSML if possible
            polyline.append(newnode)
        polylines.append(polyline)
        properties.setdefault("parent-poly", []).append(parent)

    for prop in organ.iterfind('properties'):

        if len(list(prop)) == 2:  # special fix for broken .rsml # TODO <--------- Remove after repairs...
            properties.setdefault('parent-node' , []).append(-1)

        for p in prop:  # i.e length, type, etc..
            if 'value' in p.attrib:  # no value defined
                try:
                    properties.setdefault(str(p.tag), []).append(float(p.attrib['value']))
                except:
                    properties.setdefault(str(p.tag), []).append(float(np.nan))
            else:
                try:
                    properties.setdefault(str(p.tag), []).append(float(p.text))
                except:
                    properties.setdefault(str(p.tag), []).append(float(np.nan))

    for funcs in organ.iterfind('functions'):
        for fun in funcs:
            samples = [ ]
            for sample in fun.iterfind('sample'):
                if 'value' in sample.attrib:
                    samples.append(float(sample.attrib['value']))
                else:
                    try:
                        samples.append(float(sample.text))
                    except:
                        samples.append(float(np.nan))
            functions.setdefault(str(fun.attrib['name']), []).append(samples)

    pi = len(polylines) - 1
    for elem in organ.iterfind('root'):  # and all laterals
        polylines, properties, functions = parse_rsml_(elem, polylines, properties, functions, pi)

    return polylines, properties, functions


def read_rsml(name:str) -> (list, dict, dict):
    """Parses the RSML file into:

    Args:
    name(str): file name of the rsml file

    Returns: 
    (list, dict, dict):
    (1) a (flat) list of polylines, with one polyline per root
    (2) a dictionary of properties, one per root, adds "parent_poly" holding the index of the parent root in the list of polylines
    (3) a dictionary of functions     
    """
    root = ET.parse(name).getroot()
    plant = root[1][0]
    polylines = []
    properties = { }
    functions = { }
    for elem in plant.iterfind('root'):
        (polylines, properties, functions) = parse_rsml_(elem, polylines, properties, functions, -1)

    return polylines, properties, functions


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

# DEPRICATED works only with tap root systems...
# def get_parameter(polylines :list, funcs :dict, props :dict) -> (list, list, list):
#     """ Copies radii and creation times, one value per segment """
#     radii, cts, types = [], [], []
#     for i, p in enumerate(polylines):
#         for j in range(0, len(p)):
#             if 'diameter' in funcs:
#                 radii.append(funcs['diameter'][i][j] / 2)
#             else:
#                 radii.append(0.)
#             if 'emergence_time' in funcs:
#                 cts.append(funcs['emergence_time'][i][j])
#             else:
#                 cts.append(0.)
#             if 'type' in funcs:
#                 types.append(props["type"][i])
#             else:
#                 types.append(0)
#     return radii[1:], cts[1:], types[1:]


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
        plt.plot(nodes[:, 0], nodes[:, 2], color = cmap(f(prop[i])))  # x,z plot
    plt.axis('equal')
    plt.show()


def plot_multiple_rsml(polylines:list, prop:list, times):
    """Plots the polylines in y-z axis with colors given by a root property
    TODO
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
    for i, s in enumerate(segs):
        plt.plot([nodes[s[0], 1], nodes[s[1], 1]], [nodes[s[0], 2], nodes[s[1], 2]], color = cmap(f(fun[i])))
    plt.axis('equal')
    plt.show()


if __name__ == '__main__':

    polylines, properties, functions = read_rsml("../../grids/RootSystem.rsml")
    print("Properties:")
    for key, v in properties.items():
        print("\t", key, len(properties[key]))
    print("Functions:")
    for key, v in functions.items():
        print("\t", key, len(functions[key]))
    plot_rsml(polylines, properties["type"])

    nodes, segs = get_segments(polylines, properties)
    nodes = np.array(nodes)
    segs = np.array(segs, dtype = np.int64)

#     print(functions["emergence_time"]) # TODO
#     plot_segs(nodes, segs, functions["emergence_time"])  # slow
