<img src="Logo_long_white.png" alt="drawing" width="400"/>
Try Online: [CPlantBox Web App](https://cplantbox.fz-juelich.de/)

# Introduction

CPlantBox is a functional-structural plant model that is built in a modular way that can be used at several levels of complexity. CPlantBox describes the geometry of plants by their individual organs, such as roots, stems, and leaves, which evolve over time. It can model functional aspects such as water and carbon dynamics within the plant, and provides general tools to build plant soil-interaction models. To solve partial differential equations CPlantBox can use the finite volume solver DuMu<sup>x</sup> and offers simplified Python interfaces in the repository [_dumux-rosi_](https://github.com/Plant-Root-Soil-Interactions-Modelling/dumux-rosi).   

# Installation

## Linux - with Python script
This installation method requires Ubuntu >= 20.04 and Python (>= 3.7, <3.14). For CPlantBox without _dumux-rosi_, download the Python file "installCPlantBox.py", and run it:
```bash
sudo apt-get update
sudo apt-get upgrade
[ ! -d 'cpbenv' ] && python3 -m venv cpbenv &&  source cpbenv/bin/activate ||  source cpbenv/bin/activate
wget https://raw.githubusercontent.com/Plant-Root-Soil-Interactions-Modelling/CPlantBox/master/installCPlantBox.py
python3 installCPlantBox.py
```
For CPlantBox with _dumux-rosi_, download and run the Python file "installDumuxRosi_Ubuntu.py" (the file is based on the DuMu<sup>x</sup> installation file).
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

By default, your local repository will only track the master branch. This allows for a quicker download and lower use of memory. In case you would like to switch between branches, you have two options:
1. Use the switch_branch.py script:
```bash
python3 switch_branch.py <branch_to_switch_to>
```
2. Install all the branches initially, by running
```bash
python3 installDumuxRosi_Ubuntu.py full
```
You can then change branch by doing 
```bash
git checkout <branch_to_switch_to>
```

## Linux - using a conda or Python environment

1. Clone the repository (master only) by
```bash
git clone --depth 1 -b master https://github.com/Plant-Root-Soil-Interactions-Modelling/CPlantBox.git
```
or, to download all branches:
```bash
git clone https://github.com/Plant-Root-Soil-Interactions-Modelling/CPlantBox.git
```
 Initialize the repository by
```bash
cd CPlantBox
git submodule update --init --recursive
```

2. Create the an environment

For conda use
```bash
conda env create -f environment.yml
conda activate cpb
```
For Python use
```
python3 -m venv cpb
source cpb/bin/activate
pip install -r requirements.txt
```
Finally, initialize cmake and build and install CPlantBox:
```
cmake .
make install
```
3. Test the installation by running a tutorial example, e.g.:
```bash
cd tutorial/examples/
python example1a_small.py
```

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




