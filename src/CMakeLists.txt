#
# 1. Make CPlantBox library
#

include_directories(structural)
include_directories(functional)
include_directories(visualisation)

include_directories(external/tinyxml2)
# find_package (Eigen3 3.3 REQUIRED NO_MODULE)
# link_libraries (Eigen3::Eigen)

include_directories(external/) # for PiafMunch
include_directories(${PROJECT_SOURCE_DIR}/src/external/suitsparse/include/) # for PiafMunch
include_directories(${PROJECT_SOURCE_DIR}/src/external/sundials/include/) # for PiafMunch


add_library(CPlantBox SHARED
            structural/organparameter.cpp
            structural/rootparameter.cpp
			structural/seedparameter.cpp
			structural/leafparameter.cpp
			structural/stemparameter.cpp
            structural/Organ.cpp
            structural/Root.cpp
            structural/RootDelay.cpp            
            structural/Seed.cpp
			structural/Stem.cpp
			structural/Leaf.cpp
			structural/Plant.cpp
            structural/Organism.cpp
            structural/Plant.cpp            
            structural/RootSystem.cpp
            structural/MappedOrganism.cpp
     		structural/sdf.cpp
            structural/SegmentAnalyser.cpp            
            structural/tropism.cpp            
			
            functional/XylemFlux.cpp
			functional/Perirhizal.cpp
			functional/PlantHydraulicParameters.cpp
			functional/PlantHydraulicModel.cpp			

			visualisation/Quaternion.h
			visualisation/CatmullRomSpline.h
			visualisation/PlantVisualiser.h
			visualisation/PlantVisualiser.cpp
            
			external/tinyxml2/tinyxml2.cpp		
			external/aabbcc/AABB.cc
            external/gauss_legendre/gauss_legendre.cpp
)



#
# 2. Make CPlantBox Pyhthon binding
#

# find_package(pybind11 REQUIRED) 
#set(PYBIND11_PYTHON_VERSION 3.7)
find_package( PythonInterp 3.7 REQUIRED )
find_package( PythonLibs 3.7 REQUIRED )
add_subdirectory(external/pybind11)

pybind11_add_module(plantbox SHARED 
			PyPlantBox.cpp   								
			)
target_link_libraries(plantbox PUBLIC 
			CPlantBox)
			
set_target_properties(plantbox PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/)






