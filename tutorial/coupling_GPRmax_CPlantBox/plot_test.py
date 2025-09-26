import sys; sys.path.append("../.."); sys.path.append("../../src/")
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass
import plantbox as pb
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.colors

@dataclass
class PlantParameters:
    min_b: list
    max_b: list
    cell_number: list
    param_name: str
    rh_params: str
def load_cplantbox_input(file_path: str) -> tuple:
    model_inputs = file_path.split("_")
    soil_type = "_".join(model_inputs[10:12])  # 'hydrus_loam'
    gprMax_name = file_path.split("/")[-1]
    npzfile = np.load(file_path + '.npz')
    res = model_inputs[-7]
    plant_ = model_inputs[1] 
    swc = npzfile['arr_0']  # soil water content
    # swc = np.transpose(npzfile['arr_0'], (1,0,0))  # soil water content
    print(swc.shape)
    swp = npzfile['arr_1']
    rv = npzfile['arr_2']  # root volume
    rvf = npzfile['arr_3']  # root volume fraction
    cn = npzfile['arr_4']  # cell number
    domain_x = npzfile['arr_5']  # 'target' from simulation
    domain_y = npzfile['arr_6']  # 'target' from simulation
    soil_depth = (npzfile['arr_7']/-100)
    thetas = 0.43
    
    if res == 'high':
        res_value = 1.0
    elif res == 'low':
        res_value = 3.0
    else:
        raise ValueError("Invalid resolution")

    return swc, rvf, soil_depth, thetas, domain_x, domain_y, res_value, gprMax_name, plant_

def get_plant_parameters(plant_: str, res: float, soil_depth: float) -> PlantParameters:
    if plant_ == "maize":
        min_b = [-75/2, -15/2, soil_depth] #[cm]
        max_b = [75/2, 15/2, 0.]  #[cm]
        cell_number = [int(75/res), int(15/res), int(soil_depth*-1/res)]
        param_name = "Zeamays_synMRI_modified.xml"
        rh_params = 'couvreur2012'
    elif plant_ == "wheat":
        min_b = [-15/2, -3/2, soil_depth] #[cm]
        max_b = [15/2, 3/2, 0.]  #[cm]
        cell_number = [int(15/res), int(3/res), int(soil_depth*-1/res)] 
        param_name = "wheat_Morandage.xml"
        rh_params = 'wheat_Giraud2023adapted'

    return min_b, max_b, cell_number, param_name

def get_soil_inputparameters():
    eps_w = 84  # permittivity water at 10°C
    eps_s = 4  # permittivity soil matrix
    phi = 0.43  # porosity
    eps_r = 0.8 * 84  # permittivity roots (wc root = 80%)
    cond_s = 0.015  # conductivity soil
    cond_r = cond_s  # conductivity roots
    return eps_w, eps_s, phi, eps_r, cond_s, cond_r

def calculate_eps_soil_4p(swc, rvf, phi, eps_s, eps_w, eps_r):
    eps_ss = (((1 - phi) * np.sqrt(eps_s) + swc * np.sqrt(eps_w) + (phi - swc - rvf) * 1 + rvf * np.sqrt(eps_r)) ** 2)
    return eps_ss

def plotting_cplantbox_output(swc, rvf, soil_depth, thetas, domain_x, domain_y, res, eps_ss, param_name, min_b, max_b):
    fig = plt.figure(figsize=plt.figaspect(0.5))
    """ resimulate root system """
    path = "../../modelparameter/structural/rootsystem/"
    rs = pb.MappedPlant()

    rs.readParameters(path + param_name + ".xml")
    sdf = pb.SDF_PlantBox(0.99 * (max_b[0] - min_b[0]), 0.99 * (max_b[1] - min_b[1]), max_b[2] - min_b[2])
    rs.setGeometry(sdf)
    rs.setSeed(1)
    rs.initialize()
    rs.simulate(rs_age, False)
    for i in range(0, sim_time):
        rs.simulate(i)
    ana = pb.SegmentAnalyser(rs.mappedSegments())
    segs = ana.segments
    nodes = ana.nodes
    radius = ana.getParameter("radius")

    ax1 = fig.add_subplot(1, 4, 1)
    cmap = cm.seismic
    norm = matplotlib.colors.Normalize(vmin=np.min(radius), vmax=np.max(radius))
    fc = cmap(norm(radius))
    for k, s in enumerate(segs):
        s1 = segs[k]
        n1, n2 = nodes[s1.x], nodes[s1.y]
        ax1.plot([n1.x, n2.x], [n1.z, n2.z], color = fc[k,:])
    m = cm.ScalarMappable(cmap=plt.cm.seismic)
    m.set_array([])
    plt.colorbar(m, ax = ax1, label = 'Radius (cm)')
    ax1.set_aspect('equal')
    ax1.axis('off')


    ax2 = fig.add_subplot(1, 4, 2)
    swc_ = swc[:,int(domain_y/2),:]
    colors = plt.cm.jet_r(swc_)
    norm = plt.Normalize(vmin=np.min(swc_), vmax=np.max(swc_))
    ax2.imshow(np.rot90(swc_), cmap='seismic_r', norm=norm, aspect=1) 
    m = plt.cm.ScalarMappable(cmap=plt.cm.jet_r, norm=norm)
    m.set_array([])
    plt.colorbar(m, ax=ax2, label='Soil water content (-)')
    ax2.set_title("$\\theta$ (-)")
    ax2.axis('off')

    ax3 = fig.add_subplot(1, 4, 3)
    rwc_ = rvf[:,int(domain_y/2),:]
    colors = plt.cm.jet_r(rwc_)
    norm = plt.Normalize(vmin=np.min(rwc_), vmax=np.max(rwc_))
    ax3.imshow(rwc_, cmap='seismic_r', norm=norm, aspect=1) 
    m = plt.cm.ScalarMappable(cmap=plt.cm.jet_r, norm=norm)
    m.set_array([])
    plt.colorbar(m, ax=ax3, label='Root volume fraction (-)')
    ax3.set_title("RVF")
    ax3.axis('off')

    ax4 = fig.add_subplot(1, 3, 3)
    eps_ss_ = eps_ss[:,int(domain_y/2),:]
    colors = plt.cm.jet_r(eps_ss_)
    norm = plt.Normalize(vmin=np.min(eps_ss_), vmax=np.max(eps_ss_))
    ax4.imshow(eps_ss_, cmap='jet_r', norm=norm, aspect=1) 
    m = plt.cm.ScalarMappable(cmap=plt.cm.jet_r, norm=norm)
    m.set_array([])
    plt.colorbar(m, ax=ax4, label='SWC 4-phase')
    ax4.set_title("$\\epsilon_{4phase}$ (-)")
    ax4.axis('off')

    plt.show()


if __name__ == "__main__":
    outputfolder = 'results/npz_maize_high_resolution_hydrus_loam_age70'
    out_files = os.listdir(outputfolder)
    for file in out_files:
        file_path = outputfolder + '/' + (file.split('.')[0])
        swc, rvf, soil_depth, thetas, domain_x, domain_y, res, gprMax_name, plant_ = load_cplantbox_input(file_path)
        min_b, max_b, cell_number,param_name = get_plant_parameters(plant_, res, soil_depth)
        eps_w, eps_s, phi, eps_r, cond_s, cond_r = get_soil_inputparameters()
        eps_ss = calculate_eps_soil_4p(swc, rvf, phi, eps_s, eps_w, eps_r)
        plotting_cplantbox_output(swc, rvf, soil_depth, thetas, domain_x, domain_y, res, eps_ss, param_name,min_b, max_b)