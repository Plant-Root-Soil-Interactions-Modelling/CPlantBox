import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import visualisation.vtk_plot as vp

import matplotlib.pyplot as plt
import numpy as np
import math


mycp = pb.MycorrhizalPlant()
path = "../../modelparameter/structural/rootsystem/"
name = "Anagallis_femina_Leitner_2010"
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
dispersed = True
infradius = 0
for i in range(0, len(root)):
    root[i].infradius = infradius
    if root[i].infradius != 0:
        dispersed = False
        root[i].f_inf = pb.SoilLookUpSDF(infbox, 1, 0.0, 0.1)
        print("ATTENTION: The infection radius is not 0, the infection will be local")
    root[i].dx = 0.05
    

mycp.initialize(True)

# --- Simulationseinstellungen ---
simtime = 10            # Gesamtdauer (Tage)
fpd = 250                # Schritte pro Tag
N = simtime * fpd       # Gesamtanzahl Schritte
dt = simtime / N        # Schrittweite in Tagen
time = np.linspace(0, simtime, N)
ratio = False            # Set to True for ratio plots, False for absolute values

# --- Datenspeicher initialisieren ---
observed_primary = []
observed_new_primary = []
expected_new_primary = []
expected_primary = []

observed_secondary = []
total_lengths = []
nonmycL = []

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

# Globaldaten
observed_primary.append(sum(mycp.getParameter("primaryInfection")))
observed_secondary.append(sum(mycp.getParameter("secondaryInfection")))
total_lengths.append(sum(mycp.getParameter("length")))
nonmycL.append(total_lengths[0] - observed_primary[0] - observed_secondary[0])

# Erwartungswertberechnung
P = mycp.getOrganRandomParameter(pb.root)[1].p * dt
expected_primary.append(0)
expected_segment_primary.append(0)

# --- Hauptzeitschleife ---
for t in range(1, len(time)):
    mycp.simulate(dt, False)
    infs = mycp.getNodeInfections(2)
    temp_ana = pb.SegmentAnalyser(mycp)

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

    observed_primary.append(sum(mycp.getParameter("primaryInfection")))
    observed_secondary.append(sum(mycp.getParameter("secondaryInfection")))
    total_lengths.append(sum(mycp.getParameter("length")))

    # --- Differenzen und Erwartungswerte ---
    delta_obs = observed_primary[t] - observed_primary[t - 1]
    delta_obs_segment = observed_segment_primary[t] - observed_segment_primary[t - 1]
    observed_new_primary.append(delta_obs)
    observed_segment_new_primary.append(delta_obs_segment)

    nonmycL.append(total_lengths[t] - observed_primary[t] - observed_secondary[t])
    delta_expected = nonmycL[t] * P
    delta_expected_segment = observed_segment_noninf[t] * P

    expected_new_primary.append(delta_expected)
    expected_primary.append(expected_primary[-1] + delta_expected)

    expected_segment_new_primary.append(delta_expected_segment)
    expected_segment_primary.append(expected_segment_primary[-1] + delta_expected_segment)
    # --- Ausgabe der Wahrscheinlichkeiten --- 
    obtf = math.floor(N/10)  # Beobachtungszeitraum für die empirische Wahrscheinlichkeit
    if t % obtf == 0 and t > 0:
        # Beobachtet: Gesamtlänge
        total_observed = observed_primary[t] - observed_primary[t-obtf]
        total_nonmyc = sum(nonmycL[t-obtf:t])
        # print(mycp.getOrganRandomParameter(pb.root)[1].p)
        P_empirisch = total_observed / total_nonmyc if total_nonmyc > 0 else 0

        # Beobachtet: Segmentbasiert
        total_observed_seg = observed_segment_primary[t] - observed_segment_primary[t-obtf]
        total_nonmyc_seg = sum(observed_segment_noninf[t-obtf:t])
        P_empirisch_seg = total_observed_seg / total_nonmyc_seg if total_nonmyc_seg > 0 else 0

        print(f"Zeitschritt t={t} mit Beobachtunszeitraum deltat={obtf}:")
        print(f"- Definiert:      P = {P:.5f}")
        print(f"- Gemessen (gesamt):     P_emp = {P_empirisch:.5f}")
        print(f"- Gemessen (Segment):    P_emp_seg = {P_empirisch_seg:.5f}")
        print("f totale Länge:", total_lengths[t])
        print("f total Länge Segment:", observed_segment_total[t])
        print("-" * 40)

# --- Plot ---
if ratio:
    fig, axes = plt.subplots(1, 2, figsize=(14, 5), sharex=True)

    # Maximaler Y-Wert über beide Plots hinweg
    y_max = max(
        max(observed_new_primary + expected_new_primary),
        max(observed_segment_new_primary + expected_segment_new_primary)
    ) * 1.05  # 5 % Puffer

    # Plot 1: Gesamtlängenvergleich
    axes[0].plot(time[1:], observed_new_primary, label="Beobachtet (gesamt)", color="tab:blue")
    axes[0].plot(time[1:], expected_new_primary, label="Erwartet (gesamt)", color="tab:orange")
    axes[0].set_title("Primäre Infektion (gesamt)")
    axes[0].set_xlabel("Zeit [Tage]")
    axes[0].set_ylabel("Änderung Infizierte Länge [cm]")
    axes[0].set_ylim(0, y_max)
    axes[0].legend()
    axes[0].grid(True)

    # Plot 2: Segmentbasierter Vergleich
    axes[1].plot(time[1:], observed_segment_new_primary, label="Beobachtet (Segment)", color="tab:green")
    axes[1].plot(time[1:], expected_segment_new_primary, label="Erwartet (Segment)", color="tab:red")
    axes[1].set_title("Primäre Infektion pro Segment")
    axes[1].set_xlabel("Zeit [Tage]")
    axes[1].set_ylabel("Änderung Infizierte Länge [cm]")
    axes[1].set_ylim(0, y_max)
    axes[1].legend()
    axes[1].grid(True)

    plt.tight_layout()
    plt.show()
else:
    plt.plot(np.asarray(observed_primary), label="Primary Infection")
    plt.plot(np.asarray(expected_primary), label="Expected Primary Infection")
    # plt.plot(np.asarray(observed_secondary), label="Secondary Infection")
    plt.plot(np.asarray(total_lengths), label="Length")
    # plt.plot(np.asarray(total_lengths) - np.asarray(observed_primary) - np.asarray(observed_secondary), label="Non-mycorrhizal Length")
    plt.legend()
    plt.title("Infection over time")
    plt.xlabel("Time")
    plt.ylabel("[cm]")
    plt.show()