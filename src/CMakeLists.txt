#
# 1. Make RootBox library
#

# file(GLOB all_sources *.cpp) # i changed that to exclude PythonRootSystem.cpp
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
            RootSystem.cpp
            Plant.cpp
            
            SegmentAnalyser.cpp            
            sdf.cpp
            tropism.cpp
            
			external/tinyxml2/tinyxml2.cpp
            external/aabbcc/AABB.cc)
            #external/gauss_legendre/gauss_legendre.cpp


set_target_properties(CPlantBox PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/)

#
# 2. Make plantbox library
#

# Find necessary packages for plantbox
find_package(PythonLibs 3)
find_package(PythonInterp 3)    # creates PYTHON_VERSION_MAJOR, PYTHON_VERSION_MINOR

# Only build the library if Python 3 was found
if (PythonLibs_FOUND)
  include_directories( ${PYTHON_INCLUDE_DIRS} )
#
   find_package( Boost COMPONENTS python${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR} REQUIRED )
#	find_package( Boost COMPONENTS python-py${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR} REQUIRED ) # if above is not working
  include_directories( ${Boost_INCLUDE_DIR} )

  # make shared library for Pyhton binding
  add_library(plantbox SHARED
            organparameter.cpp
            rootparameter.cpp
			leafparameter.cpp
			seedparameter.cpp
			stemparameter.cpp
            
            SegmentAnalyser.cpp
            Organ.cpp
            Root.cpp
            Seed.cpp
			Stem.cpp
			Leaf.cpp
			            
            Organism.cpp
            Plant.cpp
            RootSystem.cpp
            
            SegmentAnalyser.cpp
            sdf.cpp
            tropism.cpp

			external/tinyxml2/tinyxml2.cpp
            external/aabbcc/AABB.cc
            #external/gauss_legendre/gauss_legendre.cpp

            PyBoostPlantBox.cpp
)
  target_link_libraries(plantbox ${Boost_LIBRARIES} )
  set_target_properties(plantbox PROPERTIES PREFIX "" )
  set_target_properties(plantbox PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/)
else ()
  message(WARNING "Python 3 was not found on your system! You can only use the C++ interface of CPlantBox. \
                   Install Python 3 and rerun CMake to build the plantbox library.")
endif ()
