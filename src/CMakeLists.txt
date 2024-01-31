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
add_executable(cpb 
			cpb.cpp                        
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
            external/gauss_legendre/gauss_legendre.cpp)	
set_target_properties(CPlantBox PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/)



