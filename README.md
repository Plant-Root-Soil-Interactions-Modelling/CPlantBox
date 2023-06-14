

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
If you have ubuntu >= 20.04, you can download the python file "installDumuxRosi_Ubuntu.py".
This file is based on the dumux installation file. 
run
```bash
python3 installDumuxRosi_Ubuntu.py
```
It will create a "DUMUX" folder and install inside the dependencies necessary to run dumux-rosi.
This script might work on other linux OS but has not been tested.

If you want to install CPlantBox without Dumux, download the python file "installCPlantBox.py".
This installaiton files requires a conda environment with python >= 3.7

run
```bash
conda create -n cpb_py39 python=3.9
conda activate cpb_py39 
python3 installCPlantBox.py
```
It will create a "CPB" folder and install inside the dependencies necessary to run CPlantBox.
# Folder sructure

`/modelparameter`		Plant parameter files\
`/src`			CPlantBox C++ codes\
`/test`   Python tests for all CPlantBox classes\
`/tutorial` 		learn to use CPlantBox
`/experimental`		Specific applications (in sub-folders). contrary to scripts in `/tutorial`, might not be kept up to date\

# Code documentation

Create the documentation by running doxygen in the folder 
$ doxygen doxy_config

The documentation will be located in the folder /doc. Compile doc/latex/refman.tex to generate the full doxygen documentation in doc/latex/refman.pdf.

Collaboration diagrams give an overview of the code in folder /docs.

# Examples
Simulation videos availabe in Youtube Channel https://www.youtube.com/channel/UCPK-pFfpK94jiamgwHxX32Q

[![Plant Simulations -- 8K resolution](https://media.giphy.com/media/LmBztw7mNwluJPJ3cU/giphy.gif)](https://www.youtube.com/watch?v=jNbvjW-WFvk "CPlantBox Simulations -- 8K resolution")



