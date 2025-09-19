"""
specialized scipt for the wine rsml data 
"""
import sys; sys.path.append("../.."); sys.path.append("../../src/")
sys.path.append("../../gui/estimate/")
import visualisation.vtk_plot as vp
from estimate_data import EstimateDataModel
import estimate_plots as ep

import matplotlib.pyplot as plt
import numpy as np

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

"""
pick path and files 
"""

file_path = "../../modelparameter/structural/rootsystem"# "RSML_Fichtl_11:23"
file_names = ["B-23_Fichtl.rsml"]#["B-08.rsml", "D-25.rsml"]  # , "E-39.rsml", D-25.rsml # repetitions of same genotype

"""
RSMLs:
Type 1: the initially planted roots
Type 2: first order roots (begin to grow at once)
Type 3: second order roots (emergence times are estimated based on parametrisation of first order roots) 
"""

measurement_times = [90.] * len(file_names)  # set measurement times [days]

""" 
    open rsml files 
"""
data = EstimateDataModel()  # new data model
data.open_files(file_path, file_names)  # open rmsl
data.times = measurement_times  # initialize
for i in range(0, len(data.times)):
    data.estimates[i] = {}

""" 
    get individual root lenths from RSML files 
"""
root_lengths = [np.array(data.rsmls[i].properties["length"]) for i in range(0, len(file_names))]
indices_type2 = data.pick_order(2)
print("e.g. measurement 0 lenghts of type 2:", root_lengths[0][indices_type2[0]])

all_type2_lengths = []  # flatten all type2 lengths
for i in range(0, len(file_names)):
    all_type2_lengths.extend(root_lengths[i][indices_type2[i]])

print("Mean and SD of all Type 2 roots", np.mean(all_type2_lengths), np.std(all_type2_lengths))
# plt.hist(all_type2_lengths, 10)
# plt.xlabel("Length (cm)")
# plt.show()

""" 
    find obvious parameters 
"""
for i in range(2, 4):
    indices = data.pick_order(i)
    data.estimate_zones_(indices)  #  creates lb, ln, la, radius a, inseriton angle theta; e.g. writes single values into self.estimates[i][(j, "la")], where j is the root index
    data.aggregate_parameters_(indices, target_type = i)  # aggregates the individual root parameters (mean, sd) into data.parameters (list of RootRandomParameters) at index target_type

print()
print("****************************************")
print(data.parameters[1])  # lb, ln, la, radius a, inseriton angle theta is set
print("****************************************")

""" 
    find elongation rates for Type 3
"""
# set "age" of of order 0 to zero, and order 2 roots to measurement time
indices_type1 = data.pick_order(1)
for i, j_ in enumerate(indices_type2):
    for j in j_:
        data.estimates[i][(j, "age")] = 0.

indices_type2 = data.pick_order(2)
for i, j_ in enumerate(indices_type2):
    for j in j_:
        data.estimates[i][(j, "age")] = measurement_times[i]

# fit elongation rate for type 2
print("Type 2")
data.parameters[2].lmax = 100  # cm
data.fit_root_length_(indices_type2, base_method = 1 , target_type = 2)  # base_method = 1 fits r, base_method = 2 fits r, and lmax; note that "age" must be known for order=2 (target_type)
data.add_r_(indices_type2, 2)

# # fit elongation rate for type 3
print("\nType 3")
data.parameters[3].lmax = 100  # cm
data.compute_age(indices_type2, 2, apical_method = 1)  # sets age for the next iteration (order++), apical_method == 0, delay based; apical_method == 1, apical length based

indices_type3 = data.pick_order(3)
# for i in range(0, len(indices_type3[1])):
#     print(data.estimates[1][(indices_type3[1][i], "age")])

data.fit_root_length_(indices_type3, base_method = 1 , target_type = 3)  # base_method = 1 fits r, base_method = 2 fits r, and lmax;
# data.fit_root_length_(indices_type3, base_method = 2 , target_type = 3)

fig, ax = plt.subplots(1, 1, figsize = (9, 8))
ep.plot_laterals(data, 0, 0, ax, orders_ = [3])
plt.show()

