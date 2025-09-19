"""
Specialized scipt for the wine fine root rsml data to obtain fine root parametrs and subType=2 inter-node distances 
"""
import sys; sys.path.append("../.."); sys.path.append("../../src/")
sys.path.append("../../gui/estimate/")
import plantbox as pb

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

names_ = ["", "1st order", "2nd order", "3nd order", "higher order",
          "higher order", "higher order", "higher order", "higher order", "higher order", "higher order"]
col_ = ["y*", "r*", "b*", "g*", "c*", "m*"] * 2


def split_(values, v):  # calculates L2 on standard deviation
    std0, std1 = np.std(values[values < v]), np.std(values[values >= v])
    # print(std0 * std0 + std1 * std1)
    return std0 * std0 + std1 * std1


def split(values):
    test_values = np.linspace(np.min(values), np.max(values), 100)
    l2std = [split_(values, v) for v in test_values[1:-1]]
    split_value = test_values[np.argmin(l2std)]
    return split_value


def split_indices_ge(indices, root_lengths, v):
    indices0 = []
    for i, ind in enumerate(indices):
        ind_ = np.array(ind)
        root2 = np.array(root_lengths[i])[ind_]
        indices0.append(ind_[root2[i] >= v])
    return indices0


def plot_laterals(data, ax, indices, order):
    # plotting data
    c = np.array([len(x) for x in indices])
    print("plot_laterals(): number of order", order, "roots", c)
    age_ = []
    l_ = []
    for i, j_ in enumerate(indices):
        for j in j_:
            age_.append(data.estimates[i][(j, "age")])
            l_.append(data.rsmls[i].properties["length"][j])
    ax.plot(age_, l_, col_[order], label = names_[order])  # , alpha = i / len(data.times) * 0.8 + 0.2
    # plotting fit
    k = data.parameters[order].lmax
    r = data.parameters[order].r
    print("plot_laterals", k, r)
    max_time = np.max(age_)
    t_ = np.linspace(0, max_time, 200)
    l_ = ep.negexp_length(t_, r, k)
    ax.plot(t_, l_, "k")
    ax.set_xlabel("Estimated root age [day]")
    ax.set_ylabel("Measured root length [cm]")
    ax.legend()

"""
pick path and files 
"""
file_path = "Grapevine Estimator/finerootloss"
file_names = ["C-1E.rsml", "C-1W.rsml", "C-2E.rsml", "C-2W.rsml"]  # , "E-39.rsml", D-25.rsml # repetitions of same genotype
rootsystem_age = [123, 125, 123, 125]  # days

"""
RSMLs:
Type 1: the initially planted roots
Type 2: first order roots (begin to grow at once)
Type 3: second order roots (emergence times are estimated based on parametrisation of first order roots) 
"""

""" 
    open rsml files 
"""
data = EstimateDataModel()  # new data model
data.open_files(file_path, file_names)  # open rmsl
data.times = rootsystem_age  # initialize
for i in range(0, len(data.times)):
    data.estimates[i] = {}
print("***********\n\n\n\n")

""" 
    split order 2 roots into long and short
"""
root_lengths = [np.array(data.rsmls[i].properties["length"]) for i in range(0, len(file_names))]
indices_type2 = data.pick_order(2)
print("e.g. measurement 0 lenghts of type 2:", root_lengths[1][indices_type2[0]])

all_type2_lengths = []  # flatten all type2 lengths
for i in range(0, len(file_names)):
    all_type2_lengths.extend(root_lengths[i][indices_type2[i]])
all_type2_lengths = np.array(all_type2_lengths)

all_median = np.median(all_type2_lengths)
print("Mean and SD of all Type 2 roots", all_median, np.mean(all_type2_lengths), np.std(all_type2_lengths))

shorter = all_type2_lengths[all_type2_lengths < all_median]
print("Mean of roots < median", np.mean(shorter), np.std(shorter), "n =", len(shorter))
longer = all_type2_lengths[all_type2_lengths >= all_median]
print("Mean of roots >= median", np.mean(longer), np.std(longer), "n =", len(longer))

print("\nlets decrease mean and sd by clustering\n")

values = all_type2_lengths
split_value2 = split(values)
print("split_index", split_value2)
shorter = values[values < split_value2]
print("Mean of roots < median", np.mean(shorter), np.std(shorter), "n =", len(shorter))
longer = values[values >= split_value2]
print("Mean of roots >= median", np.mean(longer), np.std(longer), "n =", len(longer))
print()

# print()  # do it per measurement
# for i in range(0, len(indices_type2)):
#     print("Measurement", i)
#     ind_ = np.array(indices_type2[i])
#     lengths_ = root_lengths[i][ind_]
#     split_value = split(np.array(lengths_))
#     print("split_value", split_value)
#     shorter = lengths_[lengths_ < split_value]
#     print("Mean of roots < median", np.mean(shorter), np.std(shorter), "n =", len(shorter))
#     longer = lengths_[lengths_ >= split_value]
#     print("Mean of roots >= median", np.mean(longer), np.std(longer), "n =", len(longer))

# plt.hist(all_type2_lengths, 20)
# plt.xlabel("Length (cm)")
# plt.show()

"""
    find obvious parameters
"""
for i in range(1, 4):
    indices = data.pick_order(i)
    data.estimate_zones_(indices)  #  creates lb, ln, la, radius a, inseriton angle theta; e.g. writes single values into self.estimates[i][(j, "la")], where j is the root index
    data.aggregate_parameters_(indices, target_type = i)  # aggregates the individual root parameters (mean, sd) into data.parameters (list of RootRandomParameters) at index target_type

# for i in range(1, 4):
#     print("Radius", data.parameters[i].a, data.parameters[2].a_s)

# print("\n")
# print(data.parameters[2])  # lb, ln, la, radius a, inseriton angle theta is set # TODO theta should be regarding z axis (for type 2), should work if we detach it

"""
    find elongation rates for Type 2
"""
# set "age" of of order 0 to zero, and order 2 roots to measurement time
indices_type1 = data.pick_order(1)
for i, j_ in enumerate(indices_type1):
    for j in j_:
        data.estimates[i][(j, "age")] = 0.

indices_type2 = data.pick_order(2)
for i, j_ in enumerate(indices_type2):
    for j in j_:
        data.estimates[i][(j, "age")] = rootsystem_age[i]

indices_type2 = data.pick_order(2)
all_type2_lengths = []  # flatten all type2 lengths
for i in range(0, len(file_names)):
    all_type2_lengths.extend(root_lengths[i][indices_type2[i]])
all_type2_lengths = np.array(all_type2_lengths)

indices_type3 = data.pick_order(3)
all_type3_lengths = []  # flatten all type3 lengths
for i in range(0, len(file_names)):
    all_type3_lengths.extend(root_lengths[i][indices_type3[i]])
all_type3_lengths = np.array(all_type3_lengths)

# fit elongation rate for type 2
print("\nType 2: r")
data.parameters[2].lmax = 100  # cmwriteParameters
data.parameters[2].lmax = np.max(all_type2_lengths)  # cm
indices_larger = split_indices_ge(indices_type2, root_lengths, split_value2)

# data.fit_root_length_(indices_type2, base_method = 1 , target_type = 2)  # base_method = 1 fits r, base_method = 2 fits r, and lmax; note that "age" must be known for order=2 (target_type)
indices_larger = indices_type2
# print(len(indices_larger))
# print(indices_larger[0])
# print(indices_larger[1])
# print(indices_larger[2])
# print(indices_larger[3])

data.fit_root_length_(indices_larger, base_method = 1 , target_type = 2)  # base_method = 1 fits r, base_method = 2 fits r, and lmax; note that "age" must be known for order=2 (target_type)
data.add_r_(indices_type2, 2)

# fig, ax = plt.subplots(1, 1, figsize = (9, 8))
# plot_laterals(data, ax, indices = data.pick_order(2), order = 2)
# plt.show()

# fit elongation rate for type 3
print("\nType 3 r")
data.parameters[3].lmax = 100  # cm
data.parameters[3].lmax = np.max(all_type3_lengths)  # cm
data.compute_age(indices_type2, 2, apical_method = 1)  # sets age for the next iteration (order++), apical_method == 0, delay based; apical_method == 1, apical length based

indices_type3 = data.pick_order(3)
# for i in range(0, len(indices_type3[1])):
#     print(data.estimates[1][(indices_type3[1][i], "age")])

data.fit_root_length_(indices_type3, base_method = 1 , target_type = 3)  # base_method = 1 fits r, base_method = 2 fits r, and lmax;
# data.fit_root_length_(indices_type3, base_method = 2 , target_type = 3)

# fig, ax = plt.subplots(1, 1, figsize = (9, 8))
# plot_laterals(data, ax, indices = data.pick_order(3), order = 3)
# plt.show()

""" write set """
plant = pb.Plant()

for i in range(0, 4):
    data.parameters[i].subType = i
    param = data.parameters[i].copy(plant)
    # print(param.subType)
    # print(param.lmax)
    # print(param.r)
    plant.setOrganRandomParameter(param)

plant.writeParameters("wine.xml")
