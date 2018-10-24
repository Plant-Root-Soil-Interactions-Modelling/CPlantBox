# CPlantBox

Just uncomment the desired example in the `examples/main.cpp` file and compile and run it.

```bash
cmake .
make
cd examples && ./test_cplantbox
```

`cmake . ` runs CMake which configures the CRootBox libraries. `make` builds the libraries and the C++ example. `cd examples && ./test_cplantbox` runs the example.

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

To build the shared library py_rootbox for coupling with Python pleaser refer to 'python building guide.txt'

# Example

Visualized in R, thanks  [@guillaumelobet](https://github.com/guillaumelobet)

![alt text](https://github.com/Plant-Root-Soil-Interactions-Modelling/CPlantBox/blob/master/results/plant.gif "Tree with leafs")

