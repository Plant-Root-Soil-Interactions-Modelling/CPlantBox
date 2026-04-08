"""vizualisation for the manuscript figure/icons:
2D tesselation and corresponding ciruclar perirhizal zones (regular vs. random points)
"""

import random

import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import ConvexHull, Voronoi

import plantbox as pb
from plantbox.functional.Perirhizal import *


def generate_random_points(num_points, range_):
    points = []
    for _ in range(num_points):
        x = random.uniform(-range_, range_)
        y = random.uniform(-range_, range_)
        points.append(np.array([x, y]))
    return np.array(points)


def generate_regular_points(xnum_points, ynum_points, xrange, yrange):
    points = []
    for i in range(ynum_points):
        for j in range(xnum_points):
            x = (j + 0.5) * xrange / xnum_points - xrange / 2
            y = (i + 0.5) * yrange / ynum_points - yrange / 2
            points.append(np.array([x, y]))
    return np.array(points)


def flip_(nodes, center, axis):
    """flips the nodes around the center according axis"""
    # print(nodes.shape)
    # print(center.shape)
    # print(np.ones((nodes.shape[0], 1)).shape)
    n_ = nodes - center
    n_[:, axis] = -n_[:, axis]
    n_ = n_ + center
    return n_


def mirror(nodes):
    """adds mirrored nodes to the 4 sides"""
    width_ = np.array([[20, 20]])
    center_ = np.array([[0, 0]])
    nodes_surr = nodes
    fipped_n = [flip_(nodes, center_, i_) for i_ in range(0, 2)]  # flip them ...
    zeros = np.zeros((nodes.shape[0], 1))
    translate_ = np.ones((nodes.shape[0], 1)) @ width_
    trans0 = np.hstack((translate_[:, 0, np.newaxis], zeros))
    trans1 = np.hstack((zeros, translate_[:, 1, np.newaxis]))
    nodes_surr = np.vstack((nodes_surr, fipped_n[0] + trans0))  # add them
    nodes_surr = np.vstack((nodes_surr, fipped_n[0] - trans0))
    nodes_surr = np.vstack((nodes_surr, fipped_n[1] + trans1))
    nodes_surr = np.vstack((nodes_surr, fipped_n[1] - trans1))
    return nodes_surr


def vizualize(points):

    figure, axes = plt.subplots()
    axes.set_aspect(1)

    # Volumes
    vol = np.empty((nodes.shape[0]))
    vol[:] = np.nan
    for idx, reg_num in enumerate(vor.point_region):
        indices = vor.regions[reg_num]
        i_ = reg_num - 1
        if idx < nodes.shape[0]:
            if -1 in indices:  # some regions can be opened
                vol[idx] = np.nan
            else:
                vol[idx] = ConvexHull(vor.vertices[indices]).volume

    # Visualize regions
    lines = []
    reg_ = vor.regions
    ver_ = vor.vertices
    for reg in reg_:
        if -1 not in reg:
            for ind in range(-1, len(reg) - 1):
                lines.append([ver_[reg[ind]], ver_[reg[ind + 1]]])
    lines = np.array(lines)

    for i in range(0, lines.shape[0]):
        plt.plot([lines[i, 0, 0], lines[i, 1, 0]], [lines[i, 0, 1], lines[i, 1, 1]], "g:")

    # Visualize Perirhizal radii in 2D
    for i in range(0, points.shape[0]):
        a = np.sqrt(vol[i] / np.pi)
        # print(vol[i], a * a * np.pi)
        circ = plt.Circle((points[i, 0], points[i, 1]), a, fill=True, alpha=0.1, linestyle="dotted", color="red")
        axes.add_artist(circ)

    plt.plot(points[:, 0], points[:, 1], "r*")
    plt.plot(ver_[:, 0], ver_[:, 1], "b*")
    plt.xlim([-10, 10])
    plt.ylim([-10, 10])

    plt.show()


# Generate 8 random 3D points
points = generate_random_points(9, 8.0)
points = generate_regular_points(3, 3, 20, 20)
print(points)


nodes = mirror(points)
vor = Voronoi(nodes)

vizualize(points)
