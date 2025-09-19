#!/bin/python3
'''
Run a factorial simulation for sugarcane root growth
@author: Mariane Dias Macedo <marianediasmacedo@gmail.com>
'''

#%%
import sys
import numpy as np
import pandas as pd
import argparse as ap
from pathlib import Path 
import itertools
import subprocess

from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank(); max_rank = comm.Get_size()

from helper_fun import is_running_in_notebook

if is_running_in_notebook():
    sys.argv = ['run_fatorial', '-D','results/fatorial/']

parser = ap.ArgumentParser(description=__doc__, formatter_class=ap.RawDescriptionHelpFormatter)
parser.add_argument('-D','--dest_path',type=str, default = 'results/fatorial/', help='Destination path for the processed data [default: results/fatorial/]')
args = parser.parse_args()

run_script_name = 'script_test.py'

# define parameter combinations
Qps = [2., 6., 10.]  # [MPa]
Ss = [55] #[55., 75., 85.]  # [%]
Treatments = ['NT31'] #['NT31', 'CT31', 'T0'] # Treatments from SRF_parameters.csv
combinations = list(itertools.product(Treatments, Qps, Ss))
print(f'Running {len(combinations)} simulations')

#%%
# Run all combinations
df_res = []
if rank < len(combinations):
    comb = combinations[rank]
    
    Treatment = comb[0]
    Qp = comb[1]
    S = comb[2]

    Qp_name = f"{Qp*1000:.0f}kpa"
    S_name = f"{S:.0f}percent"
    sim_ID = f"{Treatment}_{Qp_name}_{S_name}"
    
    cmd = f'python3 {run_script_name} -T {Treatment} -Qp {Qp} -S {S} -D {str(args.dest_path)} -P False'
    print(f'Running: {cmd}')

    subprocess.run(cmd, shell=True, check=True)

    # read csv results
    df_res = pd.read_csv(Path(args.dest_path).joinpath(f'root_tab_{sim_ID}.csv'), sep=";", decimal=",")
    #df_res.append(df_res_comb)

    # %%
    # concatenate all results and save
    # df_res = pd.concat(df_res, ignore_index=True)
    df_res.to_csv(Path(args.dest_path).joinpath('root_tab_fatorial'+str(rank)+'.csv'), sep=";", decimal=",", index=False)
    # %%
