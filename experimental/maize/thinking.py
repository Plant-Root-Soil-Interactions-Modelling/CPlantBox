import matplotlib.pyplot as plt
import numpy as np
from maize_params import *

#
# Sheat and internode lengths
#

# at least, ensure that sl_i <= il_{i+1}+sl_{i+1}


#
# Internode elongation timings and leaf emergence
#
il = get_internode_lenght()
led = get_leaf_emergence_days()
sum_il = np.sum(il)

lmax = 197  # total stem length [cm] (without mesocotyl lb)
r = lmax / 110
print("Growth rate r:", r)
print()

scaled_il = (il / sum_il) * lmax

print(sum_il, np.sum(scaled_il))

print("length")
print(np.cumsum(scaled_il))
print()

fin = np.cumsum(scaled_il) / r
print("elongation finsihed after")
print(fin)
print()

print("leaf emerge at")
print(led)
print()


print("Leaf emerge sooner than phytomere starts elongating")
print(led[1:] < fin + 2 * np.ones(fin.shape))
print()

print("Leaf emerge sooner than phytomere ends elongating")
print(led[1:-1] < fin[1:] + 2 * np.ones(fin[1:].shape))
print(fin[1:] + 2 * np.ones(fin[1:].shape))

print()


# plt.plot(np.cumsum(scaled_il) / r, np.cumsum(scaled_il), "*")
# plt.plot(led, "o")
# plt.xlabel("Time [days]")
# plt.ylabel("Cumulative length [cm]")
# plt.show()

insert_times = [8, 12, 25, 35, 45, 50, 60, 65, 70, 75, 80]
text = ["VE", "V2", "V4", "V6", "V8", "V10", "V12", "V14", "V16–V18", "VT", "R1"]
text_dec = [
    "Internodes remain very short",
    "Initial stem axis development",
    "Early internode differentiation",
    "Start of rapid basal internode elongation",
    "Basal internodes elongate rapidly",
    "Strong elongation of basal and middle internodes",
    "Middle internodes elongate rapidly",
    "Upper internodes develop; basal internodes approach final length",
    "Upper internodes elongate",
    "Plant approaches final height; most internodes fully elongated",
    "Internode elongation essentially complete",
]


# Vegetative Stages (V)

# • VE (Emergence): The coleoptile breaks through the soil surface, typically around 100–120 Growing Degree Days (GDDs).
# • V1–V(n) (Leaf Collars): Named by the number of visible leaf collars (e.g., V1 has one, V6 has six). Lower leaves may drop off later, requiring stalk dissection from the fifth node on older plants.
# • V6 (Six-Leaf Collar): The growing point moves above the soil surface, stalk elongation ramps up, and nodal roots become the primary source of water and nutrients.
# • V12 (Twelve-Leaf Collar): Total leaf number and potential kernel/ovule capacity are set, accompanied by high nutrient demand.
# • VT (Tasseling): The final branch of the tassel is fully visible, and the plant reaches maximum height just before pollen shed.

# Reproductive Stages (R)

# • R1 (Silking): Silks emerge outside the husks and receive pollen for fertilization; this stage is highly sensitive to heat and drought stress.
# • R2 (Blister): Kernels are small and white with a clear internal fluid; moisture is roughly 85%.
# • R3 (Milk): Kernels turn yellow and contain a milky white fluid as starch accumulates; moisture is around 80%.
# • R4 (Dough): Fluid thickens to a doughy paste as starch accumulation continues; moisture is roughly 70%.
# • R5 (Dent): Kernels indent at the crown as starch hardens, with a visible "milk line" advancing toward the cob; moisture is typically around 55% or lower.
# • R6 (Physiological Maturity): Maximum dry weight is reached and a black layer forms at the kernel base, cutting off nutrient flow; kernel moisture is generally 30–35%.
