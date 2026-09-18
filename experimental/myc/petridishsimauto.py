import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import plantbox.visualisation.vtk_plot as vp
from plantbox.visualisation import figure_style
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, Normalize
import time
import matplotlib as mpl
import AMFAnalysis as amf

# I changed the c++ code, so set might not be valid anymore
good_seeds = [10, 15, 25, 30, 45, 50, 70, 105, 115, 125] # these are seeds where hyphae cross the barrier, but have not been assessed in any other way for any required behaviour

def makesimulation(seed):
    mycp = pb.MycorrhizalPlant(seed)
    path = "tomatoparameters/"
    name = "TwoHyphaePlusBAS"

    start = time.perf_counter()
    animation = True
    mycp.readParameters(path + name + ".xml", fromFile = True, verbose = True)

    ## Setting up petri dish
    diameter = 9.4
    height = 1.6
    # introduce parameters for barrier and opening
    barrier_thickness = 0.16
    barrier_height = height
    opening_length = 5.0
    opening_height = 0.2

    nRings = 15
    petri_dish, small_hyphae_dish, half_dish, rings = amf.makedishes(diameter, height, barrier_thickness, barrier_height, opening_length, opening_height, root[0].a,nRings)

    # for i, ring in enumerate(rings):
    #     vp.write_container(ring, "ring_" + str(i) + ".vtp")

    # make sure to set the seed position to 0 because of the petri dish
    seed_parameter = pb.SeedRandomParameter(mycp)
    seed_parameter.seedPos.z = -height / 2 # seed is positioned in the middle of the petri dish in the z direction
    mycp.setOrganRandomParameter(seed_parameter)

    mycp.setGeometry(half_dish)
    mycp.initialize(True)

    # set up simulation times etc.
    simtime = 5
    fps = 24
    N = fps * simtime
    dt = simtime / N

    # filename = "splitpetri_dish_parametrisation_" + str(simtime)

    filename = "petridish_nowait_seed_" + str(seed)

    print("Starting simulation with seed " + str(seed) + " and filename " + filename)
    # Start simulation
    start = time.perf_counter()
    for i in range(0, N):
        if i % 10 == 0:
            print("Step " + str(i) + " of " + str(N))
        mycp.simulate(dt,False)
        if (animation):
            amf.getMycSegmentAnalyser(mycp,filename = "animation/" + filename + "_anim" +str(i) + ".vtp", stdwrite= True)    # look at roots and container
    # vp.plot_roots_and_container(mycp,half_dish)
    # input("lengths: "+ str(sum(mycp.getParameter("length"))))


    afterroots = time.perf_counter()
    # resetting some parameters for roots
    root = mycp.getOrganRandomParameters(pb.root)
    for rp in root: 
        rp.hyphalEmergenceDensity = 4
        rp.lmbd = 0.15
        mycp.setOrganRandomParameter(rp)
    # setting up hyphal parameters

    # change geometry but only for hyphae
    mycp.changeGeometry(5, petri_dish)

    # check for percentage of colonized roots
    pCol = sum(mycp.getParameter("colonizationLength")) / sum(mycp.getParameter("length"))
    # input("Initial colonization percentage: " + str(pCol*100) + "%")
    crossed_barrier = 0
    # while pCol < 0.50:
    #     mycp.simulateColonization(dt,False)
    #     N+=1
    #     pCol = sum(mycp.getParameter("colonizationLength")) / sum(mycp.getParameter("length"))
    #     for organ in mycp.getOrgans(pb.hyphae):
    #         if organ.getParameter("subType") < 3:
    #             for node in organ.getNodes():
    #                 if node.x > -barrier_thickness/2 and node.z < opening_height-barrier_height and node.y < opening_length/2 and node.y > -opening_length/2:
    #                     break


    # vp.plot_roots_and_container(mycp,petri_dish)

    # print("Simulating hyphal growth until hyphae cross the barrier")
    # print(mycp.getSimTime(), max(mycp.getParameter("creationTime")))
    while crossed_barrier <= 3:
        N+=1    
        mycp.simulate(dt,False)
        for organ in mycp.getOrgans(pb.hyphae):
            if organ.getParameter("subType") < 3:
                for node in organ.getNodes():
                    if node.x > -barrier_thickness/2 and node.z < opening_height-barrier_height and node.y < opening_length/2 and node.y > -opening_length/2:
                        crossed_barrier += 1
        # input("Colonization percentage: " + str(pCol*100) + "%, crossed barrier: " + str(crossed_barrier))
    # look at system to see how active
    # vp.plot_roots(mycp,"active")   
    crossed_time = mycp.getSimTime()
    # print(crossed_time)
    hours_hyphae = 60 ### HIER VERÄNDERUNG DAUER HYPHEN SIMULATION
    tip_densities = list()
    print(crossed_time,mycp.getSimTime(), max(mycp.getParameter("creationTime")))
    
    for i in range(0, hours_hyphae):
        if i % 10 == 0:
            print("Simulating hyphal growth step " + str(i+1) + " of " + str(hours_hyphae))
        # mycp.simulateHyphae(dt,False)
        mycp.simulate(dt,False)
        # mycp.turnOffSidePetriDish(-barrier_thickness/2,opening_height-barrier_height,  opening_length/2, -opening_length/2)
        if animation:
            ana = amf.getMycSegmentAnalyser(mycp,filename = "animation/" + filename + "_anim" +str(i) + ".vtp", stdwrite= True)            
        # print(crossed_time, mycp.getSimTime(), max(mycp.getParameter("creationTime")))
        # raise Exception
    endsim = time.perf_counter()
    # vp.plot_roots_and_container(mycp,petri_dish)
    print("Time for simulation: ", endsim-start)

    # raise Exception
    # ana = ana = amf.getMycSegmentAnalyser(mycp)

    if not animation:
        ana = amf.getMycSegmentAnalyser(mycp, write = [],stdwrite= True,filename = "animation/" + filename, step = N)

    times = np.linspace(crossed_time, max(mycp.getParameter("creationTime"))+0.01, 100)
    print('mycp.getSimTime()',mycp.getSimTime(),'max(mycp.getParameter("creationTime"))',max(mycp.getParameter("creationTime")))
    
    # times = np.linspace(0, 28, 100)

    tip_densities = amf.getParaDistperRing("nodeTips", times, ana, rings)
    ana_densities = amf.getParaDistperRing("anastomosis", times, ana, rings)
    length_densities = amf.getParaDistperRing("length", times, ana, rings)
    lengthsSubtype = amf.getParameterOverTime("length", times, ana, np.array([1,2,3]))
    times = times - crossed_time
    return tip_densities, ana_densities, length_densities, times, lengthsSubtype

diameter = 9.4
radius = diameter / 2

nRings = 15
location = np.array([radius*np.sqrt(i/nRings) for i in range(1, nRings+1)])
simulations = []

allsims = time.perf_counter()
for i in good_seeds:
    tip_densities, ana_densities, length_densities, times, lengthsSubType = makesimulation(i)
    simulations.append({
    "tip_dens": tip_densities,
    "times": times[1:],
    "lengths": lengthsSubType,
    "ana_densities": ana_densities,
    "length_densities": length_densities
})
allsims_end = time.perf_counter()
print("Time for all simulations: ", allsims_end-allsims)

# for i in range(len(good_seeds)):
#     runner_hyphae = simulations[i]["lengths"][1] + simulations[i]["lengths"][0]
#     BAS_hyphae = simulations[i]["lengths"][2]
#     ratio = BAS_hyphae / (BAS_hyphae + runner_hyphae) if (BAS_hyphae + runner_hyphae) != 0 else 0
#     print("Ratio BAS/Runner: ", ratio)
#     print(simulations[i]["lengths"])

for i,seed in enumerate(good_seeds):
    plt.pcolormesh(
        simulations[i]["times"], 
        location, 
        simulations[i]["tip_dens"], 
        shading='auto', 
        cmap='plasma'
    )
    plt.colorbar(label="Tip Density [mm$^{-2}$]")

    plt.xlabel("Time [days]")
    plt.ylabel("Distance [cm] from centre")
    # plt.title("Radial movement of hyphal tip frequency")
    plt.savefig(f"plots/tip_densities_seed_{seed:03d}.png", dpi=300, bbox_inches="tight")
    plt.close()

    plt.pcolormesh(
            simulations[i]["times"], 
            location, 
            simulations[i]["ana_densities"], 
            shading='auto', 
            cmap='plasma'
        )
    plt.colorbar(label="Tip Density [mm$^{-2}$]")
    
    plt.xlabel("Time [days]")
    plt.ylabel("Distance [cm] from centre")
        # plt.title("Radial movement of hyphal tip frequency")
    plt.savefig(f"plots/ana_densities_seed_{seed:03d}.png", dpi=300, bbox_inches="tight")
    plt.close()



t_common = np.linspace(
    max(sim["times"][0] for sim in simulations),
    min(sim["times"][-1] for sim in simulations),
    99
)
from scipy.interpolate import interp1d

tip_dens_interp = []

for sim in simulations:
    f = interp1d(sim["times"], sim["tip_dens"], axis=1)
    tip_dens_interp.append(f(t_common))

tip_dens_interp = np.stack(tip_dens_interp)   # (n_sim, n_r, n_t)

mean = tip_dens_interp.mean(axis=0)
std  = tip_dens_interp.std(axis=0)

# mean = np.mean([sim["tip_dens"] for sim in simulations], axis=0)
# std  = np.std([sim["tip_dens"] for sim in simulations], axis=0)
# t_common = simulations[0]["times"]

indices = [i*7 for i in range(15)]

fig, ax = plt.subplots(figsize=(6,5))

cmap = plt.get_cmap("plasma")
colors = cmap(np.linspace(0.2, 0.9, len(indices)))

for color, i in zip(colors, indices):

    ax.plot(
        location,
        mean[:, i],
        color=color,
        lw=2,
        label=f"t = {t_common[i]:.2f}"
    )

    ax.fill_between(
        location,
        mean[:, i] - std[:, i],
        mean[:, i] + std[:, i],
        color=color,
        alpha=0.3
    )
ax.set_xlabel(r"Distance from origin $r$ (mm)")
ax.set_ylabel(r"Tip amount)")
# ax.set_ylabel(r"Anastomosis frequency (mm$^{-2}$)")
ax.legend()

plt.tight_layout()
plt.show()
##############################
# print(np.array(tip_densities).reshape((-1, len(tip_densities[0]))))
## The problem is that the tips should be more evenly distributed i.e. more rings should have tips in them. but right now just a handful do
# tip_densities = np.array(tip_densities).reshape((-1, len(tip_densities[0])))
# tip_densities = np.transpose(tip_densities)


# Z = np.array(tip_densities)   # (r, t)

# n_r, n_t = Z.shape
# times_arr = np.array(times[1:])  # muss 99 lang sein
# fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11,4))

# cmap = plt.get_cmap("plasma")
# norm = mpl.colors.Normalize(vmin=times_arr.min(), vmax=times_arr.max())

# sm = mpl.cm.ScalarMappable(cmap=cmap, norm=norm)
# sm.set_array([])

# for i, t in enumerate(times_arr):
#     ax2.plot(
#         location,
#         Z[:, i],
#         color=cmap(norm(t)),
#         alpha=0.8,
#         lw=1
#     )

# ax2.set_xlabel("Radius r (mm)")
# ax2.set_ylabel("Tip density (mm$^{-2}$)")
# ax2.set_title("Active tip density over time")

# threshold = 0.2 * np.max(Z)
# aligned_curves = []

# for i in range(Z.shape[1]):

#     curve = Z[:, i]

#     # arrival = first position where activity starts
#     idx = np.argmax(curve > threshold)

#     arrival_r = location[idx] if idx > 0 else location[0]

#     shifted_r = location - arrival_r

#     ax1.plot(
#         shifted_r,
#         curve,
#         color=cmap(norm(times_arr[i])),
#         alpha=0.8,
#         lw=1
#     )
# ax1.set_xlabel(r"Shifted radius $r - r_{\mathrm{arrival}}$ (mm)")
# ax1.set_ylabel("Tip density (mm$^{-2}$)")
# ax1.set_title("Wave-aligned activity")
# cbar = fig.colorbar(sm, ax=[ax1, ax2])
# cbar.set_label("Time (days)")
# # plt.tight_layout()
# endplot = time.perf_counter()
# plt.show()

# # print("Time for simulation: ", endsim-start)
# # print("Time for plotting: ", endplot-endsim)
# # print("Total time: ", endplot-start)
# # loc_grid, time_grid = np.meshgrid(location, times[1:],indexing='ij')
# # loc_flat = loc_grid.ravel()
# # time_flat = time_grid.ravel()
# # tips_flat = tip_densities.ravel()

# # print("Lengths: ", len(loc_flat), len(time_flat), len(tips_flat))

# # cmap = plt.get_cmap('viridis')

# # plt.figure(figsize=(9, 5))

# # extent = [
# #     times[1], times[-1],   # x: Zeit
# #     1, len(rings)          # y: Ringe
# # ]

# # plt.contour(location, timesuse, tip_densities, levels=50, cmap='viridis')
# # plt.xlabel("Number of Ring")
# # plt.title("Radial movement of hyphal tip density")

# # plt.show()

# # sc = plt.scatter(
# #     loc_flat,                 # x = Distanz
# #     tips_flat,                # y = Messwert (kann auch 0 sein, wenn du nur Farben willst)
# #     c=time_flat,              # Farbe = Zeitstempel
# #     cmap=cmap,
# #     s=50,                     # Marker‑Größe
# #     edgecolor='k',
# #     linewidth=0.4,
# # )

# # cbar = plt.colorbar(sc, label='Time [h]')   # Legende für die Zeitfarbe
# # plt.xlabel('distance from centre [cm]')
# # plt.ylabel('Tip Count')
# # plt.grid(True, ls='--', alpha=0.5)
# # plt.tight_layout()
# # plt.show()

# # for ti, t in enumerate(time_grid):
# #     if ti % 2 == 0:  # nur jede 10. Zeitstufe plotten, um Überladung zu vermeiden
# #         plt.plot(
# #             loc_grid,               # x‑Achse: Ort
# #             tip_densities[:, ti],                # y‑Achse: Messwerte dieser Zeit
# #             label=f't={t}s',
# #             color=cmap(ti/len(time_grid)),   # gleiche Farbskala wie beim Scatter
# #             linewidth=2,
# #             marker='o',
# #             markersize=2,
# #             # markeredgecolor='k',
# #         )

# # plt.xlabel('distance [cm]')
# # plt.ylabel('Tip Count')
# # plt.title('TipCount over distance and time')
# # plt.grid(True, ls='--', alpha=0.5)
# # plt.tight_layout()
# # plt.show()
# # vp.plot_roots_and_container(mycp,petri_dish)
# # vp.write_container(petri_dish, "petri_dish.vtp")