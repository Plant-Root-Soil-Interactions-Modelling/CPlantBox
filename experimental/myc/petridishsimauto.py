import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import plantbox.visualisation.vtk_plot as vp
import numpy as np
import matplotlib.pyplot as plt
import time
import AMFAnalysis

mycp = pb.MycorrhizalPlant(100)
path = "tomatoparameters/"
name = "TwoHyphaePlusBAS"

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


## Setting up petri dish
diameter = 9.4
radius = diameter / 2
height = 1.6
# introduce parameters for barrier and opening
barrier_thickness = 0.16
barrier_height = height
opening_length = 5.0
opening_height = 0.2


# make sure to set the seed position to 0 because of the petri dish
seed_parameter = pb.SeedRandomParameter(mycp)
seed_parameter.seedPos.z = -height / 6 # seed is positioned in the middle of the petri dish in the z direction
seed_parameter.seedPos.x = - 1.0
seed_parameter.seedPos.y = 0
mycp.setOrganRandomParameter(seed_parameter)

dishes = AMFAnalysis.PetriDishSetup(diameter, height, barrier_thickness, barrier_height, opening_length, opening_height)
petri_dish = dishes[0]
half_dish = dishes[1]
moved_helper_dish_hyphae = dishes[2]
mycp.setGeometry(half_dish)
mycp.initialize(True)

# set up simulation times etc.
simtime = 5
fps = 12
N = fps * simtime
dt = simtime / N

filename = "splitpetri_dish" + str(simtime)
# Start simulation
start = time.perf_counter()
for i in range(0, N):
    if i % 6 == 0:
        print("Step " + str(i) + " of " + str(N))
    mycp.simulate(dt,True)
    if (animation):
        ana = AMFAnalysis.getMycSegmentAnalyser(mycp)
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

hours_hyphae = 120
tip_densities = list()
for i in range(0, hours_hyphae):
    print("Simulating hyphal growth step " + str(i+1) + " of 60")
    mycp.simulateHyphae(dt,False)
    AMFAnalysis.activeCondition(mycp, barrier_thickness, opening_height, barrier_height, opening_length)
    organs = mycp.getOrgans()
    for organ in organs:
        stayactive = False
        for node in organ.getNodes():
            if node.x > -barrier_thickness/2 and node.z < opening_height-barrier_height and node.y < opening_length/2 and node.y > -opening_length/2:
                stayactive = True
        organ.setActive(stayactive)
    ana = AMFAnalysis.getMycSegmentAnalyser(mycp)
    if animation:
        ana.crop(half_dish)
        ana.write("animation/" + filename + "_" + str(N+i) + ".vtp", ["radius", "subType", "creationTime", "organType", "infection", "infectionTime", "anastomosis"])
end = time.perf_counter()
ana = AMFAnalysis.getMycSegmentAnalyser(mycp)

if not animation:
    ana.write(filename + str(N+i) + ".vtp", ["radius", "subType", "creationTime", "organType", "infection", "infectionTime", "anastomosis","nodeTips"])

# set the observation "rings"
nRings = 50
mean_x = np.mean([node.x for node in nodes_crossed])
mean_y = np.mean([node.y for node in nodes_crossed])
rings = AMFAnalysis.EquiAreaRings(nRings, radius, True, mean_x, mean_y, moved_helper_dish_hyphae, height)

times = np.linspace(0, max(mycp.getParameter("creationTime"))+0.01, 100)

tip_densities = AMFAnalysis.getParaDistperRing("nodeTips", times, ana, rings, True)
# print(np.array(tip_densities).reshape((-1, len(tip_densities[0]))))
## The problem is that the tips should be more evenly distributed i.e. more rings should have tips in them. but right now just a handful do
# tip_densities = np.array(tip_densities).reshape((-1, len(tip_densities[0])))
# tip_densities = np.transpose(tip_densities)

location = [radius*np.sqrt(i/nRings) for i in range(1, nRings+1)]

loc_grid, time_grid = np.meshgrid(location, times[1:],indexing='ij')
loc_flat = loc_grid.ravel()
time_flat = time_grid.ravel()
tips_flat = tip_densities.ravel()

print("Lengths: ", len(loc_flat), len(time_flat), len(tips_flat))

cmap = plt.get_cmap('viridis')

plt.figure(figsize=(9, 5))

sc = plt.scatter(
    loc_flat,                 # x = Distanz
    tips_flat,                # y = Messwert (kann auch 0 sein, wenn du nur Farben willst)
    c=time_flat,              # Farbe = Zeitstempel
    cmap=cmap,
    s=50,                     # Marker‑Größe
    edgecolor='k',
    linewidth=0.4,
)

cbar = plt.colorbar(sc, label='Time [h]')   # Legende für die Zeitfarbe
plt.xlabel('distance from centre [cm]')
plt.ylabel('Tip Count')
plt.grid(True, ls='--', alpha=0.5)
plt.tight_layout()
plt.show()

for ti, t in enumerate(time_grid):
    if ti % 5 == 0:  # nur jede 10. Zeitstufe plotten, um Überladung zu vermeiden
        plt.plot(
            loc_grid,               # x‑Achse: Ort
            tip_densities[:, ti],                # y‑Achse: Messwerte dieser Zeit
            label=f't={t}s',
            color=cmap(ti/len(time_grid)),   # gleiche Farbskala wie beim Scatter
            linewidth=2,
            marker='o',
            markersize=2,
            # markeredgecolor='k',
        )

plt.xlabel('Distanz / Ort [m]')
plt.ylabel('Messwert')
plt.title('Messwerte über Distanz – unterschiedliche Zeiten als farbige Kurven')
# plt.legend(title='Zeitpunkt', bbox_to_anchor=(1.02, 1), loc='upper left')
plt.grid(True, ls='--', alpha=0.5)
plt.tight_layout()
plt.show()
# vp.plot_roots_and_container(mycp,rings[99])
# vp.write_container(petri_dish, "petri_dish.vtp")