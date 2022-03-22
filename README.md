
# Hi, I am CPlantBox
[![Plant Simulations -- 8K resolution](https://media.giphy.com/media/LmBztw7mNwluJPJ3cU/giphy.gif)](https://www.youtube.com/watch?v=jNbvjW-WFvk "CPlantBox Simulations -- 8K resolution")

## For simulations in https://doi.org/10.1101/810507 please use branch isp

[![DOI](https://zenodo.org/badge/95107851.svg)](https://zenodo.org/badge/latestdoi/95107851) [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/xiaoranzhou/cpb-binder/master)
## I can :
1. Create multiple plant structure
2. Coupling with PiafMunch, and make carbon and water flow inside of the plant.

## Try me 1 click

The most convenient way is to use google colab, which is a Linux virtual machine with jupyter notebook interface.
You can click the link to follow the guide there, just to click some buttons and you will be able to create plants
[here is the link to use it](http://b.cplantbox.com).

## or build local
### windows
CPlantBox is currently not available on windows. 
Some pointers to setup a linux environment on windows are given on the [wiki](https://github.com/Plant-Root-Soil-Interactions-Modelling/CPlantBox/wiki/Help-for-windows-users).
### linux
 Run CMake which configures the CPlantBox libraries by 
```bash
cmake . && make
```
in the root folder, and run some Python tutorial examples (see tutorial/latex/PlantBox_RootSytem), e.g
```bash
cd tutorial/examples/python
python3 example1a.py
```

### Ubuntu >= 20.04
If you have ubuntu >= 20.04, you can download the python file "installDumuxRosiWithUbuntu2004.py".
This file is based on the dumux installation file. 
run
```bash
python3 installDumuxRosiWithUbuntu2004.py
```
It will create a "DUMUX" folder and install inside the dependencies necessary to run dumux-rosi.
This script might work on other linux OS but has not been tested.

# Folder sructure

`/modelparameter`		Plant parameter files\
`/src`			CPlantBox C++ codes\
`/test`   Python tests for all CPlantBox classes\
`/tutorial` 		learn to use CPlantBox

# Code documentation

Create the documentation by running doxygen in the folder 
$ doxygen doxy_config

The documentation will be located in the folder /doc. Compile doc/latex/refman.tex to generate the full doxygen documentation in doc/latex/refman.pdf.

Collaboration diagrams give an overview of the code in folder /docs.

# Examples
Simulation videos availabe in Youtube Channel https://www.youtube.com/channel/UCPK-pFfpK94jiamgwHxX32Q

[![Plant Simulations -- 8K resolution](https://media.giphy.com/media/LmBztw7mNwluJPJ3cU/giphy.gif)](https://www.youtube.com/watch?v=jNbvjW-WFvk "CPlantBox Simulations -- 8K resolution")



