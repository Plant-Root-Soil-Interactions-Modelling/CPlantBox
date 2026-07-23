# Temperature-dependent root growth modeling in winter wheat

This repository contains the analysis and visualization code accompanying:

Hack J., Ullah S., Leitner D., Vanderborght J., Kage H., Schnepf A. XXXX. Soil temperature controls seasonal root system development in virtual field simulations with CPlantBox. XXX

Please cite the above publication when using this code.


## Overview

The repository contains the following scripts:

ruthe_global_variables.py - Specifies all variables necessary for the simulation as well as paths to input files and output directories.

ruthe_simulation.py - main simulation loop for temperature dependent growth. Depends on:  
    - mini_rhizotron_simulation.py - Simulates virtual rhizotubes & extracts root scores.  
    - soil_core_simulation.py - Simulates virtual soil cores & extracts root length densities.  
    - simulation_utils.py - Temperature-scaling and utility functions used by both simulation scripts.  
    - analysis_and_plotting_utils.py - Handles plotting and calculates model evaluation metrics.  

volume_weighted_voronoi.py - Simulates a single plant in a periodic domain. The three-dimensional structure is used for Voronoi tessellation and to calculate volume-weighted perirhizal radii. The root length density (RLD) per layer and the parameters of the simulated CV-RLD relationship curve are used to estimate perirhizal radii. Finally, both perirhizal radii distributions are plotted into the same graph and the plots saved under results/simulated_perirhizal_radii.

additional_analysis_graphs.ipynb - Plots the CV-RLD relationship for the field data as well as the simulated data and calcuates the relative standard error for the field data.

temperature_graphs.ipynb - Creates various graphs of the air and soil temperature for the different years and plots the temperature scaling curve.

soil_core_geometry.py - Visualizes the position of the soil cores within the virtual field as a sanity check.

visualize_virtual_field.py - This is the script that produced figure 2 of the publication.

visualize_single_root_system.py - Creates a 3D plot of a single root system.


### Additional scripts in the spatial_resolution_experiment folder:

relative_standard_error.ipynb - Calculates the relative standard error (RSE) for multiple soil core numbers and finally plots the RSE against the number of soil cores.

resolution_cv_rld.ipynb - Plots the CV-RLD relationship for different soil core numbers and layer thicknesses.

Furthermore, in results/ there is an additional folder called spatial_resolution_tests which contains depth profiles, perirhizal radii distributions as well as CSV files with the simulated mean and standard deviation values for various combinations of soil core numbers and layer thicknesses. The minirhizotron portion of the simulation has been disabled for this analysis.


## Data structure

The data/ folder contains the input data used by the simulations and analyses:

data_raw/ conains:  
    - CSV files of the soil core and minirhizotron data from the field  
    - Excel files of the climate data for 1994-98  
    - XML files with the root system architectural parameters

data_processed/ contains:  
    - pickle files of the soil temperature for the three growing periods  
    - CSV versions of the climate data Excel files  
    - cleaned climate data CSVs


## Common simulation modifications

### Disabling temperature scaling

Temperature scaling can be disabled for both the soil core as well as the minirhizotron part of the simulation. For this, the following lines have to be commented out in either mini_rhizotron_simulation.py, soil_core_simulation.py or both.  
In the function initialize_root_system comment out these lines:  
"""  
for root_type in root_system.getRootRandomParameter():  
    root_type.f_se = scale_elongation  
"""

In the function run_daily_growth_simulation comment out these lines:  
"""  
soil_temp = initialize_soil_temperature(soil_temp_df, time_step)  
scales = compute_temperature_scaling(soil_temp, dataset)  
scale_elongation.data = scales  
"""

### Changing the number of virtual soil cores within spatial-resolution-experiment

The number of virtual soil cores sampled per simulation run is controlled by NUMBER_OF_SOIL_CORE_SAMPLES in ruthe_global_variables.py. The layer thickness is controlled by NUMBER_OF_LAYERS. Please note that this is the number of layers, so the thickness is calculated as 120cm/NUMBER_OF_LAYERS.


## External Dependencies

### Python
The project requires Python 3 and the dependencies listed in pyproject.toml.

In addition, the package python3-tk is required for certain visualization and graphical components.

### CPlantBox
This repository depends on CPlantBox, which must be installed separately (On Linux or in WSL as CPlantBox does not support native Windows).
Repository: https://github.com/Plant-Root-Soil-Interactions-Modelling/CPlantBox

For full reproducibility, please use the following CPlantBox commit:

240a0642

Using a different version of CPlantBox may lead to differences in simulation results.


## Reproducing the results

All analysis scripts and figure-generation workflows used in the publication are included in this repository.

The recommended workflow is:

## Contact

For questions regarding the code or the publication, please contact:

Julia Hack
j.hack@fz-juelich.de