import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import visualisation.vtk_plot as vp

import matplotlib.pyplot as plt
import numpy as np
import math

input("Only run this script if secondary infection is disabled.\nPress Enter to continue...")
mycp = pb.MycorrhizalPlant()
path = "../../modelparameter/structural/rootsystem/"
name = "Glycine_max"
# name = "Heliantus_Pagès_2013"
mycp.readParameters(path + name + ".xml", fromFile = True, verbose = True)

hyphae_parameter = pb.HyphaeRandomParameter(mycp)
hyphae_parameter.subType = 1
hyphae_parameter.a = 0.01
hyphae_parameter.v = 1
mycp.setOrganRandomParameter(hyphae_parameter)
# print(hyphae_parameter)

root = mycp.getOrganRandomParameter(pb.root)
for rp in root:
    rp.hyphalEmergenceDensity = 1;

infbox = pb.SDF_PlantBox(3, 3, 3)
infbox = pb.SDF_RotateTranslate(infbox, 0, 0, pb.Vector3d(0, 0, -10))
local = False
for i in range(0, len(root)):
    if local:
        root[i].f_inf = pb.SoilLookUpSDF(infbox, 1, 0.0, 0.1)
        print("ATTENTION: The infection radius is not 0, the infection will be local")
    root[i].dx = 0.05
    



# --- Simulationseinstellungen ---
simtime = 10            # Gesamtdauer (Tage)
fpd = 25                # Schritte pro Tag
N = simtime * fpd       # Gesamtanzahl Schritte
dt = simtime / N        # Schrittweite in Tagen
time = np.linspace(0, simtime, N)
ratio = False            # Set to True for ratio plots, False for absolute values

def runsimulation(seed):
    mycp.setSeed(seed)
    mycp.initialize(True)
    # --- Datenspeicher initialisieren ---
    observed_segment_primary = []
    observed_segment_new_primary = []
    expected_segment_new_primary = []
    expected_segment_primary = []

    observed_segment_secondary = []
    observed_segment_noninf = []
    observed_segment_total = []

    # --- Initialer Schritt ---
    mycp.simulate(dt, False)
    infs = mycp.getNodeInfections(2)
    temp_ana = pb.SegmentAnalyser(mycp)
    # Segmentbasierte Infektionslängen
    primary_seg = secondary_seg = noninf_seg = 0
    for i in range(1, len(infs)):
        seg_len = temp_ana.getSegmentLength(i - 1)
        if infs[i] == 1:
            primary_seg += seg_len
        elif infs[i] == 2:
            secondary_seg += seg_len
        elif infs[i] == 0:
            noninf_seg += seg_len

    observed_segment_primary.append(primary_seg)
    observed_segment_secondary.append(secondary_seg)
    observed_segment_noninf.append(noninf_seg)
    observed_segment_total.append(primary_seg + secondary_seg + noninf_seg)

    # Erwartungswertberechnung
    P = mycp.getOrganRandomParameter(pb.root)[1].lmbd * dt  # P := dt * lmbd
    expected_segment_primary.append(0)

    # --- Hauptzeitschleife ---
    for t in range(1, len(time)):
        mycp.simulate(dt, False)
        infs = mycp.getNodeInfections(2)
        temp_ana = pb.SegmentAnalyser(mycp)

        seg_len_ = []
        primary_seg = secondary_seg = noninf_seg = 0
        for i in range(1, len(infs)):
            seg_len = temp_ana.getSegmentLength(i - 1)
            seg_len_.append(seg_len)  
            if infs[i] == 1:
                primary_seg += seg_len
            elif infs[i] == 2:
                secondary_seg += seg_len
            elif infs[i] == 0:
                noninf_seg += seg_len
        
        dx = np.median(seg_len_) 
        # print("dx", dx)
        
        observed_segment_primary.append(primary_seg)
        observed_segment_secondary.append(secondary_seg)
        observed_segment_noninf.append(noninf_seg)
        observed_segment_total.append(primary_seg + secondary_seg + noninf_seg)

        # --- Differenzen und Erwartungswerte ---
        delta_obs_segment = observed_segment_primary[t] - observed_segment_primary[t - 1]
        observed_segment_new_primary.append(delta_obs_segment)
        
        
        delta_expected_segment = observed_segment_noninf[t] * P*dx

        expected_segment_new_primary.append(delta_expected_segment)
        expected_segment_primary.append(expected_segment_primary[-1] + delta_expected_segment)
    return (time, observed_segment_primary, expected_segment_primary, 
            observed_segment_new_primary, expected_segment_new_primary,
            observed_segment_secondary, observed_segment_noninf, observed_segment_total)

all_obs_primary = []
all_exp_primary = []

# for seed in range(1, 50):
#     print(f"Running simulation with seed {seed}")
#     (time, observed_segment_primary, expected_segment_primary, 
#      _, _, _, _, _) = runsimulation(seed)
#     all_obs_primary.append(observed_segment_primary)
#     all_exp_primary.append(expected_segment_primary)

# # Umwandeln in Arrays
# obs_array = np.array(all_obs_primary)
# exp_array = np.array(all_exp_primary)

# # Mittelwert und Streuung
# mean_obs = np.mean(obs_array, axis=0)
# std_obs = np.std(obs_array, axis=0)
# mean_exp = np.mean(exp_array, axis=0)
# std_exp = np.std(exp_array, axis=0)

# # Plot
# plt.fill_between(time, mean_obs - std_obs, mean_obs + std_obs, color='blue', alpha=0.2, label="Observed ± SD")
# plt.plot(time, mean_obs, label="Mean Observed", color='blue')
# plt.plot(time, mean_exp, label="Mean Expected", color='orange')
# plt.fill_between(time, mean_exp - std_exp, mean_exp + std_exp, color='orange', alpha=0.2)

# plt.xlabel("Time [days]")
# plt.ylabel("Newly Infected Length [cm]")
# plt.title("Observed vs. Expected Infection Lengths")
# plt.legend()
# plt.grid(True)
# plt.show()

#Ergebnis-Speicher vorbereiten
all_infected = []
all_noninf = []

# Mehrere Seeds durchlaufen
for seed in range(1, 20):
    print(f"Running simulation with seed {seed}")
    (_, observed_segment_primary, _, 
     _, _, observed_segment_secondary, 
     observed_segment_noninf, _) = runsimulation(seed)

    infected = np.array(observed_segment_primary) + np.array(observed_segment_secondary)
    all_infected.append(infected)
    all_noninf.append(observed_segment_noninf)

# Arrays erstellen
infected_array = np.array(all_infected)
noninf_array = np.array(all_noninf)

# Mittelwert & Standardabweichung
mean_infected = np.mean(infected_array, axis=0)
std_infected = np.std(infected_array, axis=0)
mean_noninf = np.mean(noninf_array, axis=0)
std_noninf = np.std(noninf_array, axis=0)

# Plot
plt.fill_between(time, mean_infected - std_infected, mean_infected + std_infected, color='red', alpha=0.2, label="Infected ± SD")
plt.plot(time, mean_infected, label="Mean Infected", color='red')

plt.fill_between(time, mean_noninf - std_noninf, mean_noninf + std_noninf, color='green', alpha=0.2, label="Non-infected ± SD")
plt.plot(time, mean_noninf, label="Mean Non-infected", color='green')

plt.xlabel("Time [days]")
plt.ylabel("Segment length [cm]")
plt.title("Infection Dynamics across Simulations")
plt.legend()
plt.grid(True)
plt.show()