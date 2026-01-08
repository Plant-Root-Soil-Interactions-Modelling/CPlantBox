"""
Defines font sizes and suitable figure sizes 
To be able to adjust in one spot.
Ideally, matplotlib figures should have the same text sizes throughout the book

subplots11 square + small (figsize = (6, 6)), medium (figsize = (9, 9)), large (figsize = (12, 12)) (default single row, single column) 
subplots21 rectangular 1:2  figsize = (6, 12), default: two rows, single column 
subplots12 rectangular 2:1, figsize = (12, 6), default: single row, two columns 

add more functions if needed
"""

import matplotlib.pyplot as plt

SMALL_SIZE = 16
MEDIUM_SIZE = 16
BIGGER_SIZE = 16
plt.rc("font", size = SMALL_SIZE)  # controls default text sizes
plt.rc("axes", titlesize = SMALL_SIZE)  # fontsize of the axes title
plt.rc("axes", labelsize = MEDIUM_SIZE)  # fontsize of the x and y labels
plt.rc("xtick", labelsize = SMALL_SIZE)  # fontsize of the tick labels
plt.rc("ytick", labelsize = SMALL_SIZE)  # fontsize of the tick labels
plt.rc("legend", fontsize = SMALL_SIZE)  # legend fontsize
plt.rc("figure", titlesize = BIGGER_SIZE)  # fontsize of the figure title
plt.rcParams["mathtext.default"] = "regular"

# prop_cycle = plt.rcParams["axes.prop_cycle"]
# colors = prop_cycle.by_key()["color"]


def subplots11(nrows = 1, ncols = 1, **params):  # equals medium
    """ single figure medium size """
    return plt.subplots(nrows, ncols, **params, figsize = (9, 9))


def subplots11small(nrows = 1, ncols = 1, **params):
    """ single figure small size (font will appear bigger at fixed lenght)"""
    return plt.subplots(nrows, ncols, **params, figsize = (6, 6))


def subplots11medium(nrows = 1, ncols = 1, **params):
    """ single figure medium size """
    return plt.subplots(nrows, ncols, **params, figsize = (9, 9))


def subplots11large(nrows = 1, ncols = 1, **params):
    """ single figure large size (font will appear smaller at fixed length)"""
    return plt.subplots(nrows, ncols, **params, figsize = (12, 12))


def subplots12(nrows = 1, ncols = 2, **params):
    """ single row, two columns """
    return plt.subplots(nrows, ncols, **params, figsize = (12, 6))


def subplots12large(nrows = 1, ncols = 2, **params):
    """ single row, two columns """
    return plt.subplots(nrows, ncols, **params, figsize = (24, 12))


def subplots21(nrows = 2, ncols = 1, **params):
    """ 2 rows, single column """
    return plt.subplots(nrows, ncols, **params, figsize = (6, 12))

