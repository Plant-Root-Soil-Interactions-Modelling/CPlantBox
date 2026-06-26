import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import plantbox.visualisation.vtk_plot as vp
from plantbox.visualisation import figure_style
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, Normalize
import time
import matplotlib as mpl

mycp = pb.MycorrhizalPlant(100)
path = "tomatoparameters/"
name = "TwoHyphaePlusBAS"

start = time.perf_counter()
animation = True
mycp.readParameters(path + name + ".xml", fromFile = True, verbose = True)

### initial root parameters
root = mycp.getOrganRandomParameter(pb.root)
for rp in root:
    rp.hyphalEmergenceDensity = 0
    rp.highresolution = 0.
    rp.dx = 0.1
    rp.lmbd = 0
    rp.maxAge = 100 # maximal infection age
    # rp.a = 0.01
    rp.tropismT = 1
    mycp.setOrganRandomParameter(rp)


def getMycSegmentAnalyser(plant):
    ana = pb.SegmentAnalyser(plant)
    ana.addData("infection", plant.getNodeInfections(2))
    ana.addData("infectionTime", plant.getNodeInfectionTime(2))
    ana.addData("anastomosis", plant.getAnastomosisPoints(5))
    ana.addData("nodeTips", plant.getNodeTips(5))
    return ana

## Setting up petri dish
diameter = 9.4
radius = diameter / 2
height = 1.6
# introduce parameters for barrier and opening
barrier_thickness = 0.16
barrier_height = height
opening_length = 5.0
opening_height = 0.2

# petri dish has a radius of 9.4 cm and a height of 1.6 cm
petri_dish = pb.SDF_PlantContainer(radius,radius,height,False)
# the helper dish is used to cut the petri dish in half, it has the same radius and height as the petri dish but is rotated and translated to cut the petri dish in half
helper_dish = pb.SDF_PlantContainer(radius,radius,height,True)
# moving the helper dish such that it halves the petri dish and removes a bit more to restrict roots to the correct side of barrier
moved_helper_dish = pb.SDF_RotateTranslate(helper_dish, 0, 0, pb.Vector3d(-(radius+barrier_thickness+root[0].a), 0, 0))
half_dish = pb.SDF_Intersection(petri_dish, moved_helper_dish)
# helper container for barrier
helper_staff = pb.SDF_PlantBox(barrier_thickness,barrier_height,diameter)
# helper container for opening in barrier
helper_staff2 = pb.SDF_PlantBox(barrier_thickness,opening_height,opening_length)
# have to  move the helper staff for the right position of the opening and barrier, the midpoint is the distance from the center of the helper staff to the center of the petri dish, the bottompoint is the distance from the center of the helper staff to the center of the petri dish in the y direction, and the helper staff is moved to the position of the opening and barrier
midpoint = opening_length / 2 - diameter / 2
bottompoint = opening_height / 2 - barrier_height / 2
# container for opening moved to the right position
helper_staff2 = pb.SDF_RotateTranslate(helper_staff2, 0,0,pb.Vector3d(0, bottompoint, midpoint))
# opening made in the barrier
helper_dish2 = pb.SDF_Difference(helper_staff, helper_staff2)
# barrier moved to the right position
moved_helper_dish2 = pb.SDF_RotateTranslate(helper_dish2, 90, pb.SDF_Axis.xaxis , pb.Vector3d(0, -radius, -height/2))
petri_dish = pb.SDF_Difference(petri_dish, moved_helper_dish2)

moved_helper_dish_hyphae = pb.SDF_RotateTranslate(helper_dish, 0, 0, pb.Vector3d(-(radius+barrier_thickness/2), 0, 0))
hyphae_petri_dish = pb.SDF_Difference(petri_dish, moved_helper_dish_hyphae)

# make sure to set the seed position to 0 because of the petri dish
seed_parameter = pb.SeedRandomParameter(mycp)
seed_parameter.seedPos.z = -height / 6 # seed is positioned in the middle of the petri dish in the z direction
seed_parameter.seedPos.x = - 1.0
seed_parameter.seedPos.y = 0
mycp.setOrganRandomParameter(seed_parameter)

mycp.setGeometry(half_dish)
mycp.initialize(True)

# set up simulation times etc.
simtime = 5
fps = 12
N = fps * simtime
dt = simtime / N

filename = "splitpetri_dish_parametrisation_" + str(simtime)
# Start simulation
start = time.perf_counter()
for i in range(0, N):
    if i % 6 == 0:
        print("Step " + str(i) + " of " + str(N))
    mycp.simulate(dt,True)
    if (animation):
        ana = getMycSegmentAnalyser(mycp)
        ana.write("animation/" + filename + str(i) + ".vtp", ["radius", "subType", "creationTime", "organType", "infection", "infectionTime", "anastomosis"])
# look at roots and container
# vp.plot_roots_and_container(mycp,half_dish)


afterroots = time.perf_counter()
# resetting some parameters for roots
for rp in root:
    rp.hyphalEmergenceDensity = 4
    rp.lmbd = 1/rp.dx
    mycp.setOrganRandomParameter(rp)
# setting up hyphal parameters

# change geometry but only for hyphae
mycp.changeGeometry(5, petri_dish)

# check for percentage of colonized roots
pCol = sum(mycp.getParameter("infectionLength")) / sum(mycp.getParameter("length"))
print("Initial infection percentage: " + str(pCol*100) + "%")
days = 0
while pCol < 0.50:
    mycp.simulateInfection(0.5,False)
    days +=0.5
    N+=12
    pCol = sum(mycp.getParameter("infectionLength")) / sum(mycp.getParameter("length"))
    print("Infection percentage: " + str(pCol*100) + "% after " + str(days) + " days.")
print("Infection reached 50% after " + str(days) + " days.")

# vp.plot_roots_and_container(mycp,petri_dish)

print("Simulating hyphal growth until hyphae cross the barrier")
crossed_barrier = 0
nodes_crossed = []
while crossed_barrier < 3:
    N+=1    
    mycp.simulateHyphalGrowth(dt,False)
    mycp.simulateHyphae(dt,False)
    for organ in mycp.getOrgans(pb.hyphae):
        if organ.getParameter("subType") == 1:
            for node in organ.getNodes():
                if node.x > -barrier_thickness/2 and node.z < opening_height-barrier_height and node.y < opening_length/2 and node.y > -opening_length/2:
                    crossed_barrier += 1
                    nodes_crossed.append(node)

mean_y = np.median([node.y for node in nodes_crossed])
print("Mean x and y of crossing points: ", mean_y)
small_dish = pb.SDF_PlantContainer(radius,radius,height,False)
small_hyphae_dish = pb.SDF_Difference(small_dish, moved_helper_dish_hyphae)

# inactivating those organs that are in the root part of the compartment
print("Inactivating those hyphae that are in the root part of the petri dish")
organs = mycp.getOrgans()
for organ in organs:
    stayactive = False
    for node in organ.getNodes():
        if node.x > -barrier_thickness/2 and node.z < opening_height-barrier_height and node.y < opening_length/2 and node.y > -opening_length/2:
            stayactive = True
    organ.setActive(stayactive)

# look at system to see how active
# vp.plot_roots(mycp,"active")   

middelsimHyphae = time.perf_counter()

hours_hyphae = 20
tip_densities = list()
for i in range(0, hours_hyphae):
    print("Simulating hyphal growth step " + str(i+1) + " of " + str(hours_hyphae))
    mycp.simulateHyphae(dt,False)
    organs = mycp.getOrgans()
    for organ in organs:
        stayactive = False
        for node in organ.getNodes():
            if node.x > -barrier_thickness/2 and node.z < opening_height-barrier_height and node.y < opening_length/2 and node.y > -opening_length/2:
                stayactive = True
        organ.setActive(stayactive)
    ana = getMycSegmentAnalyser(mycp)
    if animation:
        ana.crop(small_hyphae_dish)
        ana.write("animation/" + filename + "_" + str(N+i) + ".vtp", ["radius", "subType", "creationTime", "organType", "infection", "infectionTime", "anastomosis"])
endsim = time.perf_counter()
ana = getMycSegmentAnalyser(mycp)

if not animation:
    ana.write(filename + str(N+i) + ".vtp", ["radius", "subType", "creationTime", "organType", "infection", "infectionTime", "anastomosis","nodeTips"])

# set the observation "rings"
nRings = 50
centrepoint = [0, 1.5, 0]
small_dish = pb.SDF_PlantContainer(radius*np.sqrt(1/nRings),radius*np.sqrt(1/nRings),height,False)
ringone = pb.SDF_Difference(small_dish, moved_helper_dish_hyphae)
moved_ringone = pb.SDF_RotateTranslate(ringone, 0, 0, pb.Vector3d(centrepoint[0], centrepoint[1], centrepoint[2]))
rings = []
rings.append(moved_ringone)
for i in range(2, nRings+1):
    small_dish = pb.SDF_PlantContainer(radius*np.sqrt(i/nRings),radius*np.sqrt(i/nRings),height,False)
    small_hyphae_dish = pb.SDF_Difference(small_dish, moved_helper_dish_hyphae)
    old_dish = pb.SDF_Difference(pb.SDF_PlantContainer(radius*np.sqrt((i-1)/nRings),radius*np.sqrt((i-1) /nRings),height,False),moved_helper_dish_hyphae)
    small_hyphae_dish = pb.SDF_Difference(small_hyphae_dish,old_dish)
    moved_small_hyphae_dish = pb.SDF_RotateTranslate(small_hyphae_dish, 0, 0, pb.Vector3d(centrepoint[0], centrepoint[1], 0))
    rings.append(moved_small_hyphae_dish)

times = np.linspace(0, max(mycp.getParameter("creationTime"))+0.01, 100)

def getParaDistperRing(parameter, times, plant, rings):
    paradenmat = np.zeros((len(rings),len(times[1:])))
    for k, ring in enumerate(rings):
        ringana = pb.SegmentAnalyser(plant) # need to copy the whole plant for segment analyzer since cripping to one ring removes all information outside
        ringana.crop(ring)
        for j in range(len(times[1:])-1):
            ringana.filter("creationTime", 0, np.flip(np.asarray(times))[j])
            ringana.pack()
            distrib = ringana.getSummed(parameter)
            ringana.filter("creationTime",0,np.flip(np.asarray(times))[j+1])
            ringana.pack()
            summed = ringana.getSummed(parameter)
            paradenmat[k, len(times[1:])-1 -j] = np.array(distrib-summed).sum() 
        ringana.filter("creationTime",0,np.flip(np.asarray(times))[len(times[1:])])
        ringana.pack()
        summed = ringana.getSummed(parameter)
        paradenmat[k, -1] = np.array(summed).sum() 
    return paradenmat



tip_densities = getParaDistperRing("nodeTips", times, ana, rings)

# print(np.array(tip_densities).reshape((-1, len(tip_densities[0]))))
## The problem is that the tips should be more evenly distributed i.e. more rings should have tips in them. but right now just a handful do
# tip_densities = np.array(tip_densities).reshape((-1, len(tip_densities[0])))
# tip_densities = np.transpose(tip_densities)

location = np.array([radius*np.sqrt(i/nRings) for i in range(1, nRings+1)])
Z = np.array(tip_densities)   # (r, t)

n_r, n_t = Z.shape
times_arr = np.array(times[1:])  # muss 99 lang sein
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11,4))

cmap = plt.get_cmap("plasma")
norm = mpl.colors.Normalize(vmin=times_arr.min(), vmax=times_arr.max())

sm = mpl.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])

for i, t in enumerate(times_arr):
    ax2.plot(
        location,
        Z[:, i],
        color=cmap(norm(t)),
        alpha=0.8,
        lw=1
    )

ax2.set_xlabel("Radius r (mm)")
ax2.set_ylabel("Tip density (mm$^{-2}$)")
ax2.set_title("Active tip density over time")

threshold = 0.2 * np.max(Z)
aligned_curves = []

for i in range(Z.shape[1]):

    curve = Z[:, i]

    # arrival = first position where activity starts
    idx = np.argmax(curve > threshold)

    arrival_r = location[idx] if idx > 0 else location[0]

    shifted_r = location - arrival_r

    ax1.plot(
        shifted_r,
        curve,
        color=cmap(norm(times_arr[i])),
        alpha=0.8,
        lw=1
    )
ax1.set_xlabel(r"Shifted radius $r - r_{\mathrm{arrival}}$ (mm)")
ax1.set_ylabel("Tip density (mm$^{-2}$)")
ax1.set_title("Wave-aligned activity")
cbar = fig.colorbar(sm, ax=[ax1, ax2])
cbar.set_label("Time (days)")
plt.tight_layout()
endplot = time.perf_counter()
plt.show()

print("Time for simulation: ", endsim-start)
print("Time for plotting: ", endplot-endsim)
print("Total time: ", endplot-start)
# loc_grid, time_grid = np.meshgrid(location, times[1:],indexing='ij')
# loc_flat = loc_grid.ravel()
# time_flat = time_grid.ravel()
# tips_flat = tip_densities.ravel()

# print("Lengths: ", len(loc_flat), len(time_flat), len(tips_flat))

# cmap = plt.get_cmap('viridis')

# plt.figure(figsize=(9, 5))

# extent = [
#     times[1], times[-1],   # x: Zeit
#     1, len(rings)          # y: Ringe
# ]

# plt.contour(location, timesuse, tip_densities, levels=50, cmap='viridis')
# plt.xlabel("Number of Ring")
# plt.title("Radial movement of hyphal tip density")

# plt.show()

# sc = plt.scatter(
#     loc_flat,                 # x = Distanz
#     tips_flat,                # y = Messwert (kann auch 0 sein, wenn du nur Farben willst)
#     c=time_flat,              # Farbe = Zeitstempel
#     cmap=cmap,
#     s=50,                     # Marker‑Größe
#     edgecolor='k',
#     linewidth=0.4,
# )

# cbar = plt.colorbar(sc, label='Time [h]')   # Legende für die Zeitfarbe
# plt.xlabel('distance from centre [cm]')
# plt.ylabel('Tip Count')
# plt.grid(True, ls='--', alpha=0.5)
# plt.tight_layout()
# plt.show()

# for ti, t in enumerate(time_grid):
#     if ti % 2 == 0:  # nur jede 10. Zeitstufe plotten, um Überladung zu vermeiden
#         plt.plot(
#             loc_grid,               # x‑Achse: Ort
#             tip_densities[:, ti],                # y‑Achse: Messwerte dieser Zeit
#             label=f't={t}s',
#             color=cmap(ti/len(time_grid)),   # gleiche Farbskala wie beim Scatter
#             linewidth=2,
#             marker='o',
#             markersize=2,
#             # markeredgecolor='k',
#         )

# plt.xlabel('distance [cm]')
# plt.ylabel('Tip Count')
# plt.title('TipCount over distance and time')
# plt.grid(True, ls='--', alpha=0.5)
# plt.tight_layout()
# plt.show()
# vp.plot_roots_and_container(mycp,rings[99])
# vp.write_container(petri_dish, "petri_dish.vtp")