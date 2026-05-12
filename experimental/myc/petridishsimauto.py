import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import plantbox.visualisation.vtk_plot as vp
from plantbox.visualisation import figure_style
import numpy as np
import matplotlib.pyplot as plt
import time

mycp = pb.MycorrhizalPlant(100)
path = "tomatoparameters/"
name = "TwoHyphaePlusBAS"

animation = False
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
fps = 24
anim_time = simtime
N = fps * anim_time
dt = simtime / N

filename = "splitpetri_dish" + str(simtime)

# Start simulation
start = time.perf_counter()
for i in range(0, N):
    if i % 6 == 0:
        print("Step " + str(i) + " of " + str(N))
    mycp.simulate(dt,True)
    if (animation):
        ana = pb.SegmentAnalyser(mycp)
        ana.addData("infection", mycp.getNodeInfections(2))
        ana.addData("infectionTime", mycp.getNodeInfectionTime(2))
        ana.addData("anastomosis", mycp.getAnastomosisPoints(5))
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
while crossed_barrier < 3:
    N+=1    
    mycp.simulateHyphalGrowth(dt,False)
    mycp.simulateHyphae(dt,False)
    for organ in mycp.getOrgans(pb.hyphae):
        if organ.getParameter("subType") == 1:
            for node in organ.getNodes():
                if node.x > -barrier_thickness/2 and node.z < opening_height-barrier_height and node.y < opening_length/2 and node.y > -opening_length/2:
                    crossed_barrier += 1

        # ana = pb.SegmentAnalyser(mycp)
        # ana.addData("infection", mycp.getNodeInfections(2))
        # ana.addData("infectionTime", mycp.getNodeInfectionTime(2))
        # ana.addData("anastomosis", mycp.getAnastomosisPoints(5))
        # ana.write("animation/step" + str(N+1) + ".vtp", ["radius", "subType", "creationTime", "organType", "infection", "infectionTime", "anastomosis"])

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

hours_hyphae = 30
tip_densities = list()
for i in range(0, hours_hyphae):
    print("Simulating hyphal growth step " + str(i+1) + " of 60")
    mycp.simulateHyphae(dt,False)
    organs = mycp.getOrgans()
    for organ in organs:
        stayactive = False
        for node in organ.getNodes():
            if node.x > -barrier_thickness/2 and node.z < opening_height-barrier_height and node.y < opening_length/2 and node.y > -opening_length/2:
                stayactive = True
        organ.setActive(stayactive)
    if animation:
        ana = pb.SegmentAnalyser(mycp)
        ana.addData("infection", mycp.getNodeInfections(2))
        ana.addData("nodeTips",mycp.getNodeTips(-1))
        ana.addData("infectionTime", mycp.getNodeInfectionTime(2))
        ana.addData("anastomosis", mycp.getAnastomosisPoints(5))
        ana.crop(small_hyphae_dish)
        ana.write("animation/" + filename + "_" + str(N+i) + ".vtp", ["radius", "subType", "creationTime", "organType", "infection", "infectionTime", "anastomosis"])
end = time.perf_counter()
ana = pb.SegmentAnalyser(mycp)
ana.addData("infection", mycp.getNodeInfections(2))
ana.addData("infectionTime", mycp.getNodeInfectionTime(2))
ana.addData("anastomosis", mycp.getAnastomosisPoints(5))


if not animation:
    ana.write(filename + str(N+i) + ".vtp", ["radius", "subType", "creationTime", "organType", "infection", "infectionTime", "anastomosis"])

# set the observation "rings"
nRings = 20
radius = radius/2
small_dish = pb.SDF_PlantContainer(radius*np.sqrt(1/nRings),radius*np.sqrt(1/nRings),height,False)
rings = []
rings.append(pb.SDF_Difference(small_dish, moved_helper_dish_hyphae))
for i in range(1, nRings):
    small_dish = pb.SDF_PlantContainer(radius*np.sqrt(i/nRings),radius*np.sqrt(i/nRings),height,False)
    small_hyphae_dish = pb.SDF_Difference(small_dish, moved_helper_dish_hyphae)
    old_dish = pb.SDF_Difference(pb.SDF_PlantContainer(radius*np.sqrt((i-1)/nRings),radius*np.sqrt((i-1) /nRings),height,False),moved_helper_dish_hyphae)
    small_hyphae_dish = pb.SDF_Difference(small_hyphae_dish,old_dish)
    rings.append(small_hyphae_dish)

times = np.linspace(0, max(mycp.getParameter("creationTime"))+0.1, 100)
hld_matrix = np.zeros((len(rings), len(times[1:])))

def getParaDistperRing(parameter, times, plant, rings):
    paradenmat = np.zeros((len(rings),len(times[1:])))
    for k, ring in enumerate(rings):
        ringana = pb.SegmentAnalyser(plant) # need to copy the whole plant for segment analyzer since cripping to one ring removes all information outside
        ringana.crop(ring)
        for j in range(len(times[1:])):
            ringana.filter("creationTime", 0, np.flip(np.asarray(times))[j])
            ringana.pack()
            distrib = ringana.getSummed(parameter)
            paradenmat[k, len(times[1:])-1 -j] = np.array(distrib).sum() 
    return paradenmat

for k, ring in enumerate(rings):
    ringana = pb.SegmentAnalyser(mycp)
    ringana.crop(ring)
    for j in range(len(times[1:])):
        ringana.filter("creationTime", 0, np.flip(np.asarray(times))[j])
        ringana.pack()
        distrib = ringana.getSummed("length")
        hld_matrix[k, len(times[1:])-1 -j] = np.array(distrib).sum() / (np.pi * radius**2)

print("shape:", hld_matrix.shape)
print("min:", np.min(hld_matrix))
print("max:", np.max(hld_matrix))
print("mean:", np.mean(hld_matrix))
plt.figure(figsize=(8, 5))

extent = [
    times[1], times[-1],   # x: Zeit
    1, len(rings)          # y: Ringe
]

plt.imshow(
    hld_matrix,
    aspect='auto',
    origin='lower',
    extent=extent,
    cmap='viridis'
)

plt.colorbar(label="Hyphal Length Density")

plt.xlabel("Time")
plt.ylabel("Ring (center → outside)")
plt.title("Radial movement of hyphal length density")

plt.show()
# vp.plot_roots_and_container(mycp,hyphae_petri_dish)
# vp.write_container(petri_dish, "petri_dish.vtp")