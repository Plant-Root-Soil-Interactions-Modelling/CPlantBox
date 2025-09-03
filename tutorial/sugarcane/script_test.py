"""root growth limited by carbon source"""

# %%
import numpy as np
import pandas as pd
from pathlib import Path 

import helper_fun as hf

# import CPlantBox related modules
import sys; sys.path.append("../modules"); sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src")

import plantbox as pb
import visualisation.vtk_plot as vp

# %%
# Configure simulations

# Simulation steps
simtime = 300.  # [days]
dt = 300.
steps = round(simtime / dt)

# Assuming a constant carbon source to the root system (this could be replaced by a dynamic model)
limit_carbon = False
carbon_source = 500. # [g(Root)/day]

# Assuming a constant specific root length
root_SRL = 1750 # [cm(Root)/g(Root)]
# From doi:10.1016/j.fcr.2009.09.004 (Figure 7)

SRF = 1  # Root reduction factor

# %%
# Initialize model 

# Set up depth dependent elongation scaling function
scale_elongation = pb.EquidistantGrid1D(0, -50, 100)  # for root elongation from 0 cm to -50 cm, 100 nodes = 99 layers
soil_strength =np.ones(len(scale_elongation.data)-1) * SRF # rary soil strength that proportionaly reduce root growth by 10%
scale_elongation.data = soil_strength  # set proportionality factors
#print('scale_elongation',scale_elongation)
# Proportionally scale this function
#se = pb.ProportionalElongation()
#se.setBaseLookUp(scale_elongation)

# Instantiate root system for a maize plant
rs = pb.Plant(1)
#rs = pb.RootSystem()
# rs.setSeed(0)
name = "./Sugarcane2025"
rs.readParameters(name + ".xml")


# Set the scaling function and initialize
#for p in rs.getOrganRandomParameter(pb.root):
#for p in rs.getRootRandomParameter():
#    p.f_se = se
rs.initialize()

total_root_len = rs.getSummed("length")

# %%
# Simulation loop
for step in range(0, steps):

    print("\nSimulation step", step, steps)

    # Maximal total root length increment (cm/day)
    # carbon_source and root_SRL could be replaced by dynamic models
    maxinc = carbon_source * root_SRL

    # Simulate growth considering max root increment
    if limit_carbon:
        rs.simulate(dt, maxinc, se, True)
    else:
        rs.simulate(dt)

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
    #assert used_carbon <= carbon_source*tol, f"Mismatching carbon balance, used carbon is larger than carbon source ({round((tol-1)*100,3)}% tolerance)"
    if False:
        if step == 0:
            df_out = hf.asDataFrame(rs, step)
        else:
            df_out = pd.concat([df_out,hf.asDataFrame(rs, step)], ignore_index=True)

# %%
# Write outputs and plot
if not Path("results").exists():
    Path("results").mkdir()

rs.write("results/example_carbon.vtp")
#df_out.to_csv("results/example_carbon.csv", index=False)

#df_out.columns

#df_out.plot(x="RootSystem_Age_days", y="RootTipDepthMax_cm")

ana = pb.SegmentAnalyser()
ana.addSegments(rs)
vp.plot_roots(ana, "type")

# %%
