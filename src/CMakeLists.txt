#
# 1. Make CPlantBox library
#
add_library(CPlantBox SHARED                        
            organparameter.cpp
            rootparameter.cpp
			seedparameter.cpp
			leafparameter.cpp
			stemparameter.cpp
            Organ.cpp
            Root.cpp
            RootDelay.cpp            
            Seed.cpp
			Stem.cpp
			Leaf.cpp
			Plant.cpp
            Organism.cpp
            Plant.cpp            
            RootSystem.cpp
            MappedOrganism.cpp
            XylemFlux.cpp
     		sdf.cpp
            SegmentAnalyser.cpp            
            tropism.cpp            
			external/tinyxml2/tinyxml2.cpp
            external/aabbcc/AABB.cc
            external/gauss_legendre/gauss_legendre.cpp
			)

set_target_properties(CPlantBox PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/)

#
# 2. Make CPlantBox Pyhthon binding
#
# find_package(pybind11 REQUIRED) 
set(PYBIND11_PYTHON_VERSION 3.7)
find_package( PythonInterp 3.7 REQUIRED )
find_package( PythonLibs 3.7 REQUIRED )
add_subdirectory(external/pybind11)
pybind11_add_module(plantbox SHARED 
            sdf.cpp
            organparameter.cpp
            Organ.cpp
            Organism.cpp
            
            rootparameter.cpp
			seedparameter.cpp
			leafparameter.cpp
			stemparameter.cpp			
            
            Root.cpp
            RootDelay.cpp
            Seed.cpp
			Stem.cpp
			Leaf.cpp
			Plant.cpp
			            
            Plant.cpp            
            RootSystem.cpp
            MappedOrganism.cpp
			XylemFlux.cpp
           
            SegmentAnalyser.cpp            
            tropism.cpp
            
			external/tinyxml2/tinyxml2.cpp
            external/aabbcc/AABB.cc
            external/gauss_legendre/gauss_legendre.cpp
	    	
			PyPlantBox.cpp									
			)
			
set_target_properties(plantbox PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/)
