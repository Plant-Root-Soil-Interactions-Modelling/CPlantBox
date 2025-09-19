"""root growth limited by carbon source"""

# %%
import sys
import numpy as np
import pandas as pd
import argparse as ap
from pathlib import Path 

from helper_fun import Stress_Reduction_Function, asDataFrame, is_running_in_notebook

if is_running_in_notebook():
    sys.argv = ['script_test', '-T','NT31']

# import CPlantBox related modules
import sys; sys.path.append("../modules"); sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src")

import plantbox as pb
import visualisation.vtk_plot as vp

parser = ap.ArgumentParser(description=__doc__, formatter_class=ap.RawDescriptionHelpFormatter)
parser.add_argument('-Qp','--soil_pen_res', type=float, default = 6., help='Soil penetration resistance [MPa] (default: 6.)')
parser.add_argument('-S','--soil_rel_sat', type=float, default = 55., help='Soil saturation [%] (default: 55.)')
parser.add_argument('-T','--treatment', type=str, default = 'NT31', help='Treatment name (default: "NT31")')
parser.add_argument('-D','--dest_path',type=str, default = 'results/fatorial/', help='Destination path for the processed data [default: results/fatorial/]')
parser.add_argument('-P','--plot',type=str, default = 'True', help='Plot roots at the end [default: False]')
args = parser.parse_args()

# %%
# Configure simulations

Treatment = args.treatment
Qp = args.soil_pen_res  # [MPa]
S = args.soil_rel_sat  # [%]
dest_path = Path(args.dest_path)
plot_roots = args.plot.lower() in ('true', '1', 't', 'y', 'yes')

if not dest_path.exists():
    dest_path.mkdir(parents=True)

# Simulation steps
simtime = 300.  # [days]
dt = 3
steps = round(simtime / dt)

# Assuming a constant carbon source to the root system (this could be replaced by a dynamic model)
limit_carbon = True
carbon_source = 500. # [g(Root)/day]

# Assuming a constant specific root length
root_SRL = 1750 # [cm(Root)/g(Root)]
# From doi:10.1016/j.fcr.2009.09.004 (Figure 7)

# %%
# Get stress reduction factor from parameters file
parameters_df = pd.read_csv("SRF_parameters.csv", sep=";")
assert Treatment in parameters_df["Treatment"].unique(), f"Treatment {Treatment} not found in parameters file"
idx_parameter = parameters_df["Treatment"] == Treatment

A = parameters_df.loc[idx_parameter, ["A"]].values[0][0]
B = parameters_df.loc[idx_parameter, ["B"]].values[0][0]
C = parameters_df.loc[idx_parameter, ["C"]].values[0][0]
x0 = parameters_df.loc[idx_parameter, ["x0"]].values[0][0]
y0 = parameters_df.loc[idx_parameter, ["y0"]].values[0][0]

SRF = Stress_Reduction_Function(A, B, C, x0, y0, S, Qp)
SRF = max(0.,min(SRF, 1.0))
print(f'SRF = {SRF}')

# %%
# Initialize model 

# Set up depth dependent elongation scaling function
scale_elongation = pb.EquidistantGrid1D(0, -400, 100)  # for root elongation from 0 cm to -50 cm, 100 nodes = 99 layers
soil_strength =np.ones(len(scale_elongation.data)-1) * SRF # rary soil strength that proportionaly reduce root growth by 10%
scale_elongation.data = soil_strength  # set proportionality factors
#print('scale_elongation',scale_elongation)
# Proportionally scale this function
se = pb.ProportionalElongation()
se.setBaseLookUp(scale_elongation)

# Instantiate root system for a maize plant
rs = pb.Plant(1)
#rs = pb.RootSystem()
# rs.setSeed(0)
name = "./Sugarcane2025"
rs.readParameters(name + ".xml")


# Set the scaling function and initialize
for p in rs.getOrganRandomParameter(pb.root):
#for p in rs.getRootRandomParameter():
    p.f_se = se
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
            df_out = asDataFrame(rs, step)
        else:
            df_out = pd.concat([df_out,asDataFrame(rs, step)], ignore_index=True)

# %%
# Save results

Qp_name = f"{Qp*1000:.0f}kpa"
S_name = f"{S:.0f}percent"
sim_ID = f"{Treatment}_{Qp_name}_{S_name}"

rs_fn = f'root_arc_{sim_ID}.vtp'
rs.write(str(dest_path.joinpath(rs_fn)))

# get root architecture as dataframe
df_out = asDataFrame(rs, step)

df_out["Treatment"] = Treatment
df_out["Qp_MPa"] = Qp
df_out["S_percent"] = S
df_out["SRF"] = SRF

id_cols = ['Treatment', 'Qp_MPa', 'S_percent', 'SRF', 'CURRENT.DATE']
data_cols = [col for col in df_out.columns if col not in id_cols]
df_out = df_out[id_cols + data_cols]

df_fn = f"root_tab_{sim_ID}.csv"
df_out.to_csv(dest_path.joinpath(df_fn), sep=";", decimal="," ,index=False)

#%%
if plot_roots:
    ana = pb.SegmentAnalyser()
    ana.addSegments(rs)
    vp.plot_roots(ana, "type")

# %%
