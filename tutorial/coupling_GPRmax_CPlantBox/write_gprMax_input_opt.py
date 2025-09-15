"""Creating gprMax input file from CPlantBox Output"""
import math
import numpy as np
import csv
import os, os.path
import shutil
import pandas as pd
from scipy.io import savemat
from scipy.interpolate import RegularGridInterpolator
from scipy.interpolate import griddata
import sys
from scenario_static_setup import *

def calculate_eps_soil_4p(swc, rvf, phi, eps_s, eps_w, eps_r):
    """Calculate the 4-phase CRIM permittivity"""
    eps_ss = (((1 - phi) * np.sqrt(eps_s) + swc * np.sqrt(eps_w) + (phi - swc - rvf) * 1 + rvf * np.sqrt(eps_r)) ** 2)
    return eps_ss

def calculate_eps_soil_3p(swc, phi, eps_s, eps_w):
    """Calculate the 3-phase CRIM permittivity"""
    eps_sw = ((1 - phi) * np.sqrt(eps_s) + swc * np.sqrt(eps_w) + (phi - swc)) ** 2
    return eps_sw

def load_cplantbox_input(file):
    """Load CPlantBox input"""
    model_inputs = file.split("_")
    file_name_t = file.split("/")
    file_name = file_name_t[-1]
    grid_res = model_inputs[-8]
    soil_type = "_".join(model_inputs[10:12])  # 'hydrus_loam'
    gprMax_name_temp = file.split("/")
    gprMax_name = gprMax_name_temp[-1]
    npzfile = np.load(file+ '.npz')
    print(npzfile['arr_0'])
    swc = np.transpose(npzfile['arr_0'], (1, 2, 0))  # soil water content
    swc = np.round(swc, 2)
    swp = np.transpose(npzfile['arr_1'], (1, 2, 0))
    rv = np.transpose(npzfile['arr_2'], (1, 2, 0))  # root volume
    rvf = np.transpose(npzfile['arr_3'], (1, 2, 0))  # root volume fraction
    cn = npzfile['arr_4']  # cell number
    domain_x = npzfile['arr_5']  # 'target' from simulation
    domain_y = npzfile['arr_6']  # 'target' from simulation
    #domain_x = 1.75
    #domain_y = 1.0
    soil_height = (npzfile['arr_7']/-100)
    thetas = 0.43
    indices_rvf = np.where(rvf > 0)  # find indices of cells where roots are present
    return swc, rvf, soil_height, thetas, indices_rvf, domain_x, domain_y, grid_res, gprMax_name, file_name, cn

def get_gprmax_settings():
    """Get gprMax settings"""
    domain = 2  # when 1 use domain set below || when = 2 then use domain from CPlantBox
    c_f_Mhz = 200  # Center frequency in [MHz]
    tw = 3e-8  # time window
    PML = 0.5  # buffer in [m]
    #bcells = abc_buffer * 1000  # 500 # Number of boundary cells, 10 for absorption, 20 as buffer
    h_air = 0.5  # height air layer
    #tx_depths = [0.2, 0.4, 0.6, 0.8, 1.0, 1.2]  # depths of Tx
    tx_pos = 0  # Tx position
    rx_pos = 0.75  # Rx position
    return domain, c_f_Mhz, tw, PML, h_air, tx_pos, rx_pos
# return domain, c_f_Mhz, tw, PML, bcells, h_air, tx_depths, tx_pos, rx_pos

def get_soil_inputparameters():
    """Get soil input parameters"""
    eps_w = 84  # permittivity water at 10°C
    eps_s = 4  # permittivity soil matrix
    phi = thetas  # porosity
    eps_r = 0.8 * 84  # permittivity roots (wc root = 80%)
    cond_s = 0.015  # conductivity soil
    cond_r = cond_s  # conductivity roots
    return eps_w, eps_s, phi, eps_r, cond_s, cond_r

def write_gprmax_input_file(file, swc, rvf, soil_height, thetas, indices_rvf, domain_x, domain_y, grid_res, gprMax_name,file_name, cn):
    """Write gprMax input file"""
    domain, c_f_Mhz, tw, PML, h_air, tx_pos, rx_pos = get_gprmax_settings()
    # domain, c_f_Mhz, tw, PML, bcells, h_air, tx_depths, tx_pos, rx_pos = get_gprmax_settings()
    eps_w, eps_s, phi, eps_r, cond_s, cond_r = get_soil_inputparameters()
    # Calculate parameters for input file
    name_sw_file = 'SW_wavefile.txt'
    pulse_name = 'pulse1'
    c_f = c_f_Mhz * 10 ** 6  # Umrechnung von MHz in Hz
    d_x = domain_x / 100  # domain size in x-direction
    d_y = domain_y / 100  # domain size in y-direction
    d_z = soil_height + h_air  # domain size in z-direction
    rx_steps = 0.2
    tx_pos = 0 # Tx position
    rx_pos = 0.75 # Rx position 
    if grid_res == 'high':
        grid_size = 0.01
    else:
        grid_size = 0.03
    # Calculate soil system permittivity (4-phase equation)
    eps_ss = calculate_eps_soil_4p(swc, rvf, phi, eps_s, eps_w, eps_r)
    eps_ss = np.round(eps_ss, 2)
    # Calculate mean soil water content per layer
    swc_mean_per_layer = np.mean(np.mean(swc, axis=0), axis=0)
    # Calculate 3-phase CRIM permittivity
    eps_sw = calculate_eps_soil_3p(swc_mean_per_layer, phi, eps_s, eps_w)
    # Create a numpy array with unique values to set material properties
    eps_ss_uniq = np.unique(eps_ss)
    # Create a dictionary for root soil properties
    root_soil = {}
    for k in range(len(eps_ss_uniq)):
        key = f'root_soil_{k+1}'
        root_soil[key] = [eps_ss_uniq[k], cond_r, 1, 0]
    # Write gprMax input file
    in_filename = 'gprMax_inputs/gprmax_' + gprMax_name + '.in'
    with open(in_filename, 'w') as fid:
        fid.write(f"#title: {file_name} \n")
        fid.write(f"#domain: {d_x:4.3f} {d_y:4.3f} {d_z:4.3f}\n")
        fid.write(f"#dx_dy_dz: {grid_size:2.4f} {grid_size:2.4f} {grid_size:2.4f}\n")
        fid.write(f"#time_window: {tw:.9f}\n\n")
        for key, value in root_soil.items():
            fid.write(f"#material: {value[0]:3.3f} {value[1]:3.3f} {value[2]:3.3f} {value[3]:3.3f} {key}\n")
        fid.write("\n")
        #fid.write(f"#excitation_file: {name_sw_file} slinear extrapolate\n")
        fid.write(f"#waveform: ricker 1 {c_f:10f} my_ricker\n")
        #fid.write(f"#hertzian_dipole: z {tx_x:4.3f} {soil_height:4.3f} {rx_pos:4.3f} {pulse_name} \n")
        fid.write(f"#hertzian_dipole: z {tx_pos+PML:4.3f} {domain_x/2:4.3f}  {h_air+rx_steps:4.3f} my_ricker\n")
        fid.write(f"#rx: {tx_pos+PML+rx_pos:4.3f} {domain_x/2:4.3f} {h_air+0.2:4.3f}\n")
        fid.write(f"#src_steps: 0 0 {rx_steps}  \n")
        fid.write(f"#rx_steps: 0 0 {rx_steps} \n")
        fid.write("\n")
        fid.write(f"#box: 0 0 0 {domain_x:4.3f} {domain_y:4.3f} {h_air:4.3f} half_space  n\n")
        # Add Boxes for SWC layers
        root_x, root_y, root_z = eps_ss.shape
        start_x, start_y, start_z = 0, 0, 0
        for z in range(root_z):
            for y in range(root_y):
                for x in range(root_x):
                    # Find index where the value in eps_bulk_root matches
                    match = np.where(eps_ss_uniq == eps_ss[x, y, z])[0]
                    if match.size > 0:
                        idx = int(match[0]) + 1
                        root_soil_field = f'root_soil_{idx}'
                        x0 = start_x + x*grid_size
                        y0 = start_y + y*grid_size
                        z0 = start_z + z*grid_size
                        fid.write(f"#box: {x0:4.3f} {y0:4.3f} {z0+h_air:4.3f} {x0+grid_size:4.3f} {y0+grid_size:4.3f} {z0+grid_size+h_air:4.3f} {root_soil_field}  n\n")
        fid.write("\n")
    print('gprMax input file gprmax_' + file + '.in is written')

if __name__ == "__main__":
    outputfolder = 'results/npz_maize_high_resolution_hydrus_loam_age70_evap'
    out_files = os.listdir(outputfolder)
    for file in out_files:
        file_path = outputfolder + '/' + (file.split('.')[0])
        swc, rvf, soil_height, thetas, indices_rvf, domain_x, domain_y, grid_res, gprMax_name, file_name, cn = load_cplantbox_input(file_path)
        write_gprmax_input_file(file_path, swc, rvf, soil_height, thetas, indices_rvf, domain_x, domain_y, grid_res, gprMax_name, file_name, cn)