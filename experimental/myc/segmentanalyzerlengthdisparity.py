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
mycp.readParameters(path + name + ".xml", fromFile=True, verbose=True)

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


# --- Datenspeicher initialisieren ---
observed_primary = []
expected_primary = []

observed_secondary = []
observed_total = []
observed_not_infected = []

observed_segment_primary = []
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
observed_total.append(sum(mycp.getParameter("length")))
observed_not_infected.append(observed_total[0] - observed_primary[0] - observed_secondary[0])

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
    observed_total.append(sum(mycp.getParameter("length")))

    # --- Differenzen und Erwartungswerte ---
    delta_obs = observed_primary[t] - observed_primary[t - 1]
    delta_obs_segment = observed_segment_primary[t] - observed_segment_primary[t - 1]

    observed_not_infected.append(observed_total[t] - observed_primary[t] - observed_secondary[t])
    delta_expected = observed_not_infected[t] * P
    delta_expected_segment = observed_segment_noninf[t] * P

    expected_primary.append(expected_primary[-1] + delta_expected)

    expected_segment_primary.append(expected_segment_primary[-1] + delta_expected_segment)
    # --- Ausgabe der Wahrscheinlichkeiten --- 
    obtf = math.floor(N/10)  # Beobachtungszeitraum für die empirische Wahrscheinlichkeit
    if t % obtf == 0 and t > 0:
        # Beobachtet: Gesamtlänge
        total_observed = observed_primary[t] - observed_primary[t-obtf]
        total_nonmyc = sum(observed_not_infected[t-obtf:t])

        # Beobachtet: Segmentbasiert
        total_observed_seg = observed_segment_primary[t] - observed_segment_primary[t-obtf]
        total_nonmyc_seg = sum(observed_segment_noninf[t-obtf:t])
    
