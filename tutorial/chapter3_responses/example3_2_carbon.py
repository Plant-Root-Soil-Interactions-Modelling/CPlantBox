"""root growth limited by carbon source"""

# %%
import numpy as np

# import CPlantBox related modules
import sys; sys.path.append("../.."); sys.path.append("../../src/")
import plantbox as pb
import visualisation.vtk_plot as vp

# %%
# Configure simulations

# Simulation steps
simtime = 50  # [days]
dt = 1
steps = round(simtime / dt)

# Assuming a constant carbon source to the root system (this could be replaced by a dynamic model)
carbon_source = 30. # [g(Root)/day]

# Assuming a constant specific root length
root_SRL = 2.5 # [cm(Root)/g(Root)]

# %%
# Initialize model

# Set up depth dependent elongation scaling function
scale_elongation = pb.EquidistantGrid1D(0, -50, 100)  # for root elongation from 0 cm to -50 cm, 100 nodes = 99 layers
soil_strength =np.ones(len(scale_elongation.data)-1) * 0.9 # arbitrary soil strength that proportionaly reduce root growth by 10%
scale_elongation.data = soil_strength  # set proportionality factors

# Proportionally scale this function
se = pb.ProportionalElongation()
se.setBaseLookUp(scale_elongation)

# Instantiate root system for a maize plant
#rs = pb.Plant()
rs = pb.RootSystem()
rs.setSeed(0)
name = "../../modelparameter/structural/rootsystem/Zea_mays_4_Leitner_2014"
rs.readParameters(name + ".xml")

# Set the scaling function and initialize
#for p in rs.getOrganRandomParameter(pb.root):
for p in rs.getRootRandomParameter():
    p.f_se = se
rs.initialize()

total_root_len = rs.getSummed("length")

# %%
# Simulation loop
for step in range(0, steps):

    print("\nSimulation step", step)

    # Maximal total root length increment (cm/day)
    # carbon_source and root_SRL could be replaced by dynamic models
    maxinc = carbon_source * root_SRL

    # Simulate growth considering max root increment
    rs.simulate(dt, maxinc, se, False)

    # Root growth and carbon balance
    total_root_len_ = rs.getSummed("length")
    root_len_increment = total_root_len_ - total_root_len
    total_root_len = total_root_len_

    used_carbon = root_len_increment / root_SRL
    unnused_carbon = carbon_source - used_carbon

    print("Max roots elongation is:", round(maxinc, 1), " [cm]")
    print("Roots elongated:", round(root_len_increment, 1), " [cm]")
    print("Carbon provided:", round(carbon_source, 1), " [g]")
    print("Carbon used for growth:", round(used_carbon, 1), " [g]")

    # Carbon balance check
    tol = 1.01 # 1% tolerance
    assert used_carbon <= carbon_source*tol, f"Mismatching carbon balance, used carbon is larger than carbon source ({round((tol-1)*100,3)}% tolerance)"

# %%
# Write outputs and plot
rs.write("results/example_carbon.vtp")

ana = pb.SegmentAnalyser()
ana.addSegments(rs)
vp.plot_roots(ana, "type")

# %%
