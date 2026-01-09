"""root growth limited by carbon source"""

import numpy as np

import plantbox as pb
import plantbox.visualisation.vtk_plot as vp

# Configure simulation
sim_time = 50  # days
dt = 1  # days
n_steps = round(sim_time / dt)
carbon_source = 30.0  # Assuming a constant carbon source to the root system (g_Root/day) |\label{l3_1_carbon:SetStart}|
root_SRL = 2.5  # Assuming a constant specific root length (cm_Root)/g_Root) |\label{l3_1_carbon:SetEnd}|

# Set up depth dependent elongation scaling function
scale_elongation = pb.EquidistantGrid1D(0, -50, 100)  # |\label{l3_1_carbon:GridStart}|
soil_strength = np.ones(len(scale_elongation.data) - 1) * 0.9  #
scale_elongation.data = soil_strength  # |\label{l3_1_carbon:GridEnd}|
se = pb.ProportionalElongation()  # Proportionally scale this function
se.setBaseLookUp(scale_elongation)

# Instantiate root system for a maize plant
plant = pb.Plant()
plant.setSeed(0)
filename = "../../modelparameter/structural/rootsystem/Zea_mays_4_Leitner_2014"
plant.readParameters(filename + ".xml")

# Set the scaling function and initialize
for p in plant.getOrganRandomParameter(pb.root):  # |\label{l3_1_carbon:SefStart}|
    p.f_se = se
plant.initialize()  # |\label{l3_1_carbon:SefEnd}|

total_root_len = plant.getSummed("length")

# Simulation loop
for step in range(0, n_steps):  # |\label{l3_1_carbon:LoopStart}|
    print("\nSimulation step", step)

    # Maximal total root length increment (cm/day)
    # carbon_source, root_SRL = model() |\label{l3_1_carbon:EntryModel}|
    maxinc = carbon_source * root_SRL

    # Simulate growth considering max root increment
    plant.simulate(dt, maxinc, se, True)

    # Root growth and carbon balance
    total_root_len_ = plant.getSummed("length")
    root_len_increment = total_root_len_ - total_root_len
    total_root_len = total_root_len_

    used_carbon = root_len_increment / root_SRL
    unnused_carbon = carbon_source - used_carbon

    print("Max roots elongation is:", round(maxinc, 1), " (cm)")
    print("Roots elongated:", round(root_len_increment, 1), " (cm)")
    print("Carbon provided:", round(carbon_source, 1), " (g)")
    print("Carbon used for growth:", round(used_carbon, 1), " (g)")

    # Carbon balance check
    tol = 1.01  # 1% tolerance
    assert used_carbon <= carbon_source * tol, f"Mismatching carbon balance, used carbon is larger than carbon source ({round((tol - 1) * 100, 3)}% tolerance)"  # |\label{l3_1_carbon:LoopEnd}|

# Write outputs and plot
plant.write("results/example3_1_carbon.vtp")  # |\label{l3_1_carbon:WriteStart}|

ana = pb.SegmentAnalyser()
ana.addSegments(plant)
vp.plot_roots(ana, "subType")  # press g to save the jpg |\label{l3_1_carbon:WriteEnd}|
