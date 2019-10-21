# CPlantBox
[![DOI](https://zenodo.org/badge/95107851.svg)](https://zenodo.org/badge/latestdoi/95107851) [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Plant-Root-Soil-Interactions-Modelling/CPlantBox/blob/master/tutorial/jupyter/CPlantBox_PiafMunch_Tutorial_(include_installation).ipynb]\
Two main functions
1. Create plant structure
2. Coupling with PiafMunch to have carbon and water flow inside the structure.

# Create plant structure in 5 minutes

The most convenient way is to use google colab, which is a Linux with jupyter notebook interface.
You can click the link to follow the guide there, just to click some buttons and you will be able to create plants
[here is the link to use it](https://colab.research.google.com/github/Plant-Root-Soil-Interactions-Modelling/CPlantBox/blob/master/tutorial/jupyter/CPlantBox_PiafMunch_Tutorial_(include_installation).ipynb)


# Folder sructure

/src			CPlantBox C++ codes
/examples 		Some examples how to use the CplantBox
/modelparameter		Some root parameter, and a plant parameter files
/scripts 		Pyhthon scripts for visualization with Paraview, and Matlab scripts for parameter export
/results 		Nice result images

# Documentation

Create the documentation by running doxygen in the folder 
$ doxygen doxy_config

The documentation should now be located in the folder /doc


# Example

Visualized in R, thanks  [@guillaumelobet](https://github.com/guillaumelobet)

![alt text](https://github.com/Plant-Root-Soil-Interactions-Modelling/CPlantBox/blob/master/results/plant.gif "Tree with leafs")

