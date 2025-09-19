"""root growth limited by carbon source"""

#%% Import libs
import numpy as np

# import CPlantBox related modules
import sys; sys.path.append("../.."); sys.path.append("../../src/")
import plantbox as pb
import visualisation.vtk_plot as vp

#%% Configure simulations

# Simulation steps
simtime = 50  # [days]
dt = 1
steps = round(simtime / dt)

# Assuming a constant carbon source to the root system
carbon_source = 30. # [g(Root)/day] |\label{l3_1_carbon:SetStart}|

# Assuming a constant specific root length
root_SRL = 2.5 # [cm(Root)/g(Root)] |\label{l3_1_carbon:SetEnd}|

#%% Initialize model

# Set up depth dependent elongation scaling function
scale_elongation = pb.EquidistantGrid1D(0, -50, 100)  # |\label{l3_1_carbon:GridStart}|
soil_strength =np.ones(len(scale_elongation.data)-1) * 0.9 #
scale_elongation.data = soil_strength  # |\label{l3_1_carbon:GridEnd}|

# Proportionally scale this function
se = pb.ProportionalElongation()
se.setBaseLookUp(scale_elongation)

# Instantiate root system for a maize plant
rs = pb.Plant()
rs.setSeed(0)
name = "../../modelparameter/structural/rootsystem/Zea_mays_4_Leitner_2014"
rs.readParameters(name + ".xml")

# Set the scaling function and initialize
for p in rs.getOrganRandomParameter(pb.root): # |\label{l3_1_carbon:SefStart}|
    p.f_se = se
rs.initialize() # |\label{l3_1_carbon:SefEnd}|

total_root_len = rs.getSummed("length")

#%% Simulation loop 
for step in range(0, steps): # |\label{l3_1_carbon:LoopStart}|

    print("\nSimulation step", step)

    # Maximal total root length increment (cm/day)
    # carbon_source, root_SRL = model() |\label{l3_1_carbon:EntryModel}|
    maxinc = carbon_source * root_SRL

    # Simulate growth considering max root increment
    rs.simulate(dt, maxinc, se, True)

    # Root growth and carbon balance
    total_root_len_ = rs.getSummed("length")
    root_len_increment = total_root_len_ - total_root_len
    total_root_len = total_root_len_

    used_carbon = root_len_increment / root_SRL
    unnused_carbon = carbon_source - used_carbon

    print("Max roots elongation is:",round(maxinc,1)," [cm]")
    print("Roots elongated:",round(root_len_increment,1)," [cm]")
    print("Carbon provided:",round(carbon_source,1)," [g]")
    print("Carbon used for growth:", round(used_carbon,1)," [g]")

    # Carbon balance check
    tol = 1.01 # 1% tolerance
    assert used_carbon <= carbon_source*tol, f"Mismatching carbon balance, used carbon is larger than carbon source ({round((tol-1)*100,3)}% tolerance)" # |\label{l3_1_carbon:LoopEnd}|

#%% Write outputs and plot
rs.write("results/example_carbon.vtp") # |\label{l3_1_carbon:WriteStart}|

ana = pb.SegmentAnalyser()
ana.addSegments(rs)
vp.plot_roots(ana, "type") # press g to save the jpg |\label{l3_1_carbon:WriteEnd}|
