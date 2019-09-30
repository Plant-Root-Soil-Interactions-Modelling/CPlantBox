#
# 1. Make RootBox library
#
add_library(CPlantBox SHARED                        
            organparameter.cpp
            rootparameter.cpp
			seedparameter.cpp
			leafparameter.cpp
			stemparameter.cpp
            Organ.cpp
            Root.cpp
            Seed.cpp
			Stem.cpp
			Leaf.cpp
			Plant.cpp
            Organism.cpp
            Plant.cpp            
            RootSystem.cpp
     		sdf.cpp
            SegmentAnalyser.cpp            
            tropism.cpp            
			external/tinyxml2/tinyxml2.cpp
            external/aabbcc/AABB.cc
            # external/gauss_legendre/gauss_legendre.cpp
			)

set_target_properties(CPlantBox PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/)

#
# 2. Make plantbox library
#
# find_package(pybind11 REQUIRED) 
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
            Seed.cpp
			Stem.cpp
			Leaf.cpp
			Plant.cpp
			            
            Plant.cpp            
            RootSystem.cpp
           
            SegmentAnalyser.cpp            
            tropism.cpp
            
			external/tinyxml2/tinyxml2.cpp
            external/aabbcc/AABB.cc
			
			PyPlantBox.cpp									
			)
			
set_target_properties(plantbox PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/)
