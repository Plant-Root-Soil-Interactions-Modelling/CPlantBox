
# Hi, I am CPlantBox
[![Plant Simulations -- 8K resolution](https://media.giphy.com/media/LmBztw7mNwluJPJ3cU/giphy.gif)](https://www.youtube.com/watch?v=jNbvjW-WFvk "CPlantBox Simulations -- 8K resolution")

## For simulations in https://doi.org/10.1101/810507 please use branch isp

[![DOI](https://zenodo.org/badge/95107851.svg)](https://zenodo.org/badge/latestdoi/95107851) [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Plant-Root-Soil-Interactions-Modelling/CPlantBox/blob/master/tutorial/jupyter/CPlantBox_PiafMunch_Tutorial_(include_installation).ipynb) Teaching Material LBRAI2219:[![Teaching Material LBRAI2219](https://badgen.net/badge/Launch/on%20Google%20Colab/blue?icon=terminal)](https://colab.research.google.com/github/Plant-Root-Soil-Interactions-Modelling/CPlantBox/blob/master/tutorial/jupyter/CPlantBox_Lesson.ipynb)

## I can :
1. Create multiple plant structure
2. Coupling with PiafMunch, and make carbon and water flow inside of the plant.

## Try me 1 click

The most convenient way is to use google colab, which is a Linux virtual machine with jupyter notebook interface.
You can click the link to follow the guide there, just to click some buttons and you will be able to create plants
[here is the link to use it](https://colab.research.google.com/github/Plant-Root-Soil-Interactions-Modelling/CPlantBox/blob/master/tutorial/jupyter/CPlantBox_PiafMunch_Tutorial_(include_installation).ipynb).

## or build local

 Run CMake which configures the CPlantBox libraries by 
```bash
cmake . && make
```
in the root folder, and run some Python tutorial examples (see tutorial/latex/PlantBox_RootSytem), e.g
```bash
cd tutorial/examples/python
python3 example1a.py
```

# Folder sructure

`/modelparameter`		Plant parameter files\
`/src`			CPlantBox C++ codes\
`/test`   Python tests for all CPlantBox classes\
`/tutorial` 		learn to use CPlantBox

# Code documentation

Create the documentation by running doxygen in the folder 
$ doxygen doxy_config

The documentation will be located in the folder /doc

# Examples
Simulation videos availabe in Youtube Channel https://www.youtube.com/channel/UCPK-pFfpK94jiamgwHxX32Q

[![Plant Simulations -- 8K resolution](https://media.giphy.com/media/LmBztw7mNwluJPJ3cU/giphy.gif)](https://www.youtube.com/watch?v=jNbvjW-WFvk "CPlantBox Simulations -- 8K resolution")



