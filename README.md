<img src="Logo_long_white.png" alt="drawing" width="400"/>

# Introduction

CPlantBox is a functional-structural plant model that is built in a modular way that can be used at several levels of complexity. CPlantBox describes the geometry of plants by their individual organs, such as roots, stems, and leaves, which evolve over time. It can model functional aspects such as water and carbon dynamics within the plant, and it provides tools to build plant soil-interaction models. 

# Installation
## Python script (recommended)
### Linux
This installation method requires ubuntu >= 20.04 and python >= 3.7.\
For CPlantBox <ins>__without__</ins> the dumux-rosi extension, download the python file "installCPlantBox.py".\
Run
```bash
sudo apt-get update
sudo apt-get upgrade
[ ! -d 'cpbenv' ] && python3 -m venv cpbenv &&  source cpbenv/bin/activate ||  source cpbenv/bin/activate
wget https://raw.githubusercontent.com/Plant-Root-Soil-Interactions-Modelling/CPlantBox/master/installCPlantBox.py
python3 installCPlantBox.py
```
It will create a "CPB" folder and install inside the dependencies necessary to run CPlantBox.\
For CPlantBox <ins>__with__</ins> the dumux-rosi extension, download the python file "installDumuxRosi_Ubuntu.py" (based on the dumux installation file).\
run
```bash
sudo apt-get update
sudo apt-get upgrade
[ ! -d 'cpbenv' ] && python3 -m venv cpbenv &&  source cpbenv/bin/activate ||  source cpbenv/bin/activate
wget https://raw.githubusercontent.com/Plant-Root-Soil-Interactions-Modelling/CPlantBox/master/installDumuxRosi_Ubuntu.py
python3 installDumuxRosi_Ubuntu.py
```
This will create a "DUMUX" folder and install inside the dependencies necessary to run dumux-rosi.
CPlantBox is setup within the virtual environment 'cpbenv'. \
**Do not forget to reactivate the 'cpbenv' environment when using CPlantBox:**
```bash
source cpbenv/bin/activate
```
This script might work on other linux OS but has not been tested.

### windows
CPlantBox is currently not available on windows. 
Some pointers to setup a linux environment on windows are given on the [wiki](https://github.com/Plant-Root-Soil-Interactions-Modelling/CPlantBox/wiki/Help-for-windows-users).

## Manual Linux installation 
Clone the repository by running:
```bash
git clone --depth 1 -b master https://github.com/Plant-Root-Soil-Interactions-Modelling/CPlantBox.git
```
Run CMake which configures the CPlantBox libraries by 
```bash
cmake . && make
```
in the root folder, and run some Python tutorial examples (see tutorial/latex/PlantBox_RootSytem), e.g
```bash
cd tutorial/examples/python
python3 example1a.py
```

The dependecies are listed in the requirements.txt file.
## Installation on the JSC agrocluster
Refer to the wiki:\
https://github.com/Plant-Root-Soil-Interactions-Modelling/CPlantBox/wiki/CPlantBox-on-the-J%C3%BClich-Supercomputer-cluster
# Folder sructure

`/modelparameter`		Plant parameter files\
`/src`			CPlantBox C++ codes\
`/test`   Python tests for all CPlantBox classes\
`/tutorial` 		learn to use CPlantBox\
`/experimental`		Specific applications (in sub-folders). contrary to scripts in `/tutorial`, might not be kept up to date

# Code documentation

Create the documentation by running doxygen in the folder 
$ doxygen doxy_config

The documentation will be located in the folder /doc. Compile doc/latex/refman.tex to generate the full doxygen documentation in doc/latex/refman.pdf.

Collaboration diagrams give an overview of the code in folder /docs.

# Examples

## WebApps demonstration 

[this web-based application](http://cplantbox.com) designed for conduct and visualize plant growth simulations. It is part of Xiaoran Zhou's PhD thesis. the source code is avialable at <a href="https://github.com/xiaoranzhou/cpb">github-xiaoranzhou 

## Cloud base jupyter notebooks
1. [Structure definition and analysis](https://mybinder.org/v2/gh/Plant-Root-Soil-Interactions-Modelling/CPlantBox/workshop_1111?labpath=tutorial%2Fjupyter%2Fworkshop_11_11_2024%2F1_cplantbox.ipynb)
2. [water flow in CPlantBox](https://mybinder.org/v2/gh/Plant-Root-Soil-Interactions-Modelling/CPlantBox/workshop_1111?labpath=tutorial%2Fjupyter%2Fworkshop_11_11_2024%2F2_water_flux.ipynb)

## Videos
Simulation videos availabe in Youtube Channel https://www.youtube.com/channel/UCPK-pFfpK94jiamgwHxX32Q




