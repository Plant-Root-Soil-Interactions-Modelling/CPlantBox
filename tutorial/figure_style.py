"""
Defines font sizes and suitable figure sizes; 
ideally, matplotlib figures should have the same text sizes throughout the book   

2:1, 
1:2, 
4:3, 
3:4, 
1:1 (square)
1:1 (small)
1:1 (med)
1:1 (large)

"""
import matplotlib.pyplot as plt

SMALL_SIZE = 16
MEDIUM_SIZE = 16
BIGGER_SIZE = 16
plt.rc('font', size = SMALL_SIZE)  # controls default text sizes
plt.rc('axes', titlesize = SMALL_SIZE)  # fontsize of the axes title
plt.rc('axes', labelsize = MEDIUM_SIZE)  # fontsize of the x and y labels
plt.rc('xtick', labelsize = SMALL_SIZE)  # fontsize of the tick labels
plt.rc('ytick', labelsize = SMALL_SIZE)  # fontsize of the tick labels
plt.rc('legend', fontsize = SMALL_SIZE)  # legend fontsize
plt.rc('figure', titlesize = BIGGER_SIZE)  # fontsize of the figure title
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']


def subplots21(nrows = 1, ncols = 1, **params):
    return plt.subplots(nrows, ncols, **params, figsize = (12, 6))


def subplots12(nrows = 1, ncols = 1, **params):
    return plt.subplots(nrows, ncols, **params, figsize = (6, 12))


def subplots43(nrows = 1, ncols = 1, **params):
    return plt.subplots(nrows, ncols, **params, figsize = (12, 9))


def subplots34(nrows = 1, ncols = 1, **params):
    return plt.subplots(nrows, ncols, **params, figsize = (9, 12))


def subplots11(nrows = 1, ncols = 1, **params):
    return plt.subplots(nrows, ncols, **params, figsize = (9, 9))


def subplots11small(nrows = 1, ncols = 1, **params):
    return plt.subplots(nrows, ncols, **params, figsize = (6, 6))


def subplots11medium(nrows = 1, ncols = 1, **params):
    return plt.subplots(nrows, ncols, **params, figsize = (9, 9))


def subplots11large(nrows = 1, ncols = 1, **params):
    return plt.subplots(nrows, ncols, **params, figsize = (12, 12))

