<img src="Logo_long_white.png" alt="drawing" width="400"/>

# Introduction

CPlantBox is a functional-structural plant model that is built in a modular way that can be used at several levels of complexity. CPlantBox describes the geometry of plants by their individual organs, such as roots, stems, and leaves, which evolve over time. It can model functional aspects such as water and carbon dynamics within the plant, and provides gerneral tools to build plant soil-interaction models. To solve partial differential equations CPlantBox can use the finite volume solver DuMu<sup>x</sup> and offers simplified Python interfaces in the repository _dumux-rosi_.   

# Installation

## Linux - with Python script
This installation method requires Ubuntu >= 20.04 and Python >= 3.7. For CPlantBox without _dumux-rosi_, download the Python file "installCPlantBox.py", and run it:
```bash
sudo apt-get update
sudo apt-get upgrade
[ ! -d 'cpbenv' ] && python3 -m venv cpbenv &&  source cpbenv/bin/activate ||  source cpbenv/bin/activate
wget https://raw.githubusercontent.com/Plant-Root-Soil-Interactions-Modelling/CPlantBox/master/installCPlantBox.py
python3 installCPlantBox.py
```
For CPlantBox with _dumux-rosi_, download and run the Python file "installDumuxRosi_Ubuntu.py" (the file is based on the DuMu$^x$ installation file).
```bash
sudo apt-get update
sudo apt-get upgrade
[ ! -d 'cpbenv' ] && python3 -m venv cpbenv &&  source cpbenv/bin/activate ||  source cpbenv/bin/activate
wget https://raw.githubusercontent.com/Plant-Root-Soil-Interactions-Modelling/CPlantBox/master/installDumuxRosi_Ubuntu.py
python3 installDumuxRosi_Ubuntu.py
```
The script will install DuMu<sup>x</sup> and CPlantBox, and CPlantBox is setup within the virtual environment 'cpbenv'. 
Activate the 'cpbenv' environment when using CPlantBox:
```bash
source cpbenv/bin/activate
```
The scripts might work on other Linux OS but has not been tested.

## Linux - with conda environment

This installation method uses ```conda``` to setup the building environment for CPlantBox. It'll pull the packages from the ```conda-forge``` channel to avoid licensing restrictions from default channels. For more info on conda restrictions, check this [article](https://www.fz-juelich.de/en/rse/the_latest/the-anaconda-is-squeezing-us)

1. Clone the repository:

```bash
git clone --depth 1 -b master https://github.com/Plant-Root-Soil-Interactions-Modelling/CPlantBox.git
```

2. Create the conda environment and build CPlantBox:

```bash
cd CPlantBox
conda env create -f environment.yml
conda activate cpb
git submodule update --init --recursive
cmake .
make
```

3. Test the installation by running a tutorial example, e.g.:

```bash
cd tutorial/examples/
python example1a_small.py
```

## Linux - manual installation 
Clone the repository by running
```bash
git clone --depth 1 -b master https://github.com/Plant-Root-Soil-Interactions-Modelling/CPlantBox.git
```
and use CMake to configure and compile the CPlantBox libraries 
```bash
cmake . && make
```
To test the installation run a tutorial example, e.g
```bash
cd tutorial/examples/python
python3 example1a.py
```
Dependecies are listed in the requirements.txt file.

## Windows
CPlantBox is currently not available on windows. Some pointers to setup a Linux environment on windows are given on the [wiki](https://github.com/Plant-Root-Soil-Interactions-Modelling/CPlantBox/wiki/Help-for-windows-users).

## Installation on the JSC agrocluster
Please refer to the [wiki](https://github.com/Plant-Root-Soil-Interactions-Modelling/CPlantBox/wiki/CPlantBox-on-the-J%C3%BClich-Supercomputer-cluster)

# Folder sructure

`/modelparameter`		Plant parameter files\
`/src`			CPlantBox C++ codes\
`/test`   Python tests for all CPlantBox classes\
`/tutorial` 		learn to use CPlantBox\
`/experimental`		Specific applications (in sub-folders). contrary to scripts in `/tutorial`, might not be kept up to date

# Code documentation

Create the documentation by running doxygen in the docs/ folder 
```bash
doxygen doxy_config
```
The documentation will be created in this folder. Compile doc/latex/refman.tex to generate the full doxygen documentation in doc/latex/refman.pdf. Additionally, collaboration diagrams give an overview of the code in folder /docs.

# Online resources

## WebApps

The official [CPlantBox webapp](https://cplantbox.fz-juelich.de) helps to demonstrate the impact of various CPlantBox parameters and to analyse and explore the resulting 3D plant geometry.  

Another [web application](http://cplantbox.com) was designed to conduct simulations and visualize the dynamics of plant growth. The source code is avialable at [github-xiaoranzhou](https://github.com/xiaoranzhou/cpb).

## Jupyter notebooks
Direct links to notebooks from the summer school of 2025:
1. [Structure definition and analysis](https://mybinder.org/v2/gh/Plant-Root-Soil-Interactions-Modelling/CPlantBox/ss2025?urlpath=%2Fdoc%2Ftree%2Ftutorial%2Fjupyter%2Fsummer_school_2025%2F1_root_architecture.ipynb)
2. [Water flow in CPlantBox](https://mybinder.org/v2/gh/Plant-Root-Soil-Interactions-Modelling/CPlantBox/ss2025?urlpath=%2Fdoc%2Ftree%2Ftutorial%2Fjupyter%2Fsummer_school_2025%2F2_root_hydraulics.ipynb)
3. [Online repository](https://mybinder.org/v2/gh/Plant-Root-Soil-Interactions-Modelling/CPlantBox/ss2025_binder)

For the third repository, you need to go the folder `dumux/CPlantBox/tutorial/jupyter/summer_school_2025` to access the notebooks.

## Videos
Simulation videos availabe in [Youtube channel](https://www.youtube.com/channel/UCPK-pFfpK94jiamgwHxX32Q).




