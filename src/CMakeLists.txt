#
# 1. Make PlantBox library
#

# file(GLOB all_sources *.cpp) # i changed that to exclude PythonRootSystem.cpp
add_library(CPlantBox
              analysis.cpp
              ModelParameter.cpp
              Plant.cpp
			  Organ.cpp
			  Root.cpp
              sdf.cpp
              RootTropism.cpp
			  StemTropism.cpp
			  LeafTropism.cpp
              Stem.cpp
			  Leaf.cpp
			  Seed.cpp
			  tinyxml2.cpp
            #../examples/Exudation/gauss_legendre.cpp
			)

set_target_properties(CPlantBox PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)

#
# 2. Make py_rootbox library
#

# Find necessary packages for py_plantbox
find_package(PythonLibs 3)
find_package(PythonInterp 3) # creates PYTHON_VERSION_MAJOR, PYTHON_VERSION_MINOR

# Only build the library if Python 3 was found
if (PythonLibs_FOUND)
  include_directories( ${PYTHON_INCLUDE_DIRS} )
  find_package( Boost COMPONENTS python-py${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR} REQUIRED )
  include_directories( ${Boost_INCLUDE_DIR} )

  # make shared library for Pyhton binding
  add_library(py_plantbox SHARED
              PythonPlant.cpp
              analysis.cpp
              ModelParameter.cpp
              Plant.cpp
			  Organ.cpp
			  Root.cpp
              sdf.cpp
              RootTropism.cpp
			  StemTropism.cpp
			  LeafTropism.cpp
              Stem.cpp
			  Leaf.cpp
			  Seed.cpp
			  tinyxml2.cpp
			  
			  #../examples/Exudation/gauss_legendre.cpp
			  )
  target_link_libraries(py_plantbox ${Boost_LIBRARIES} )
  set_target_properties(py_plantbox PROPERTIES PREFIX "" )
  set_target_properties(py_plantbox PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/python)
else ()
  message(WARNING "Python 3 was not found on your system! You can only use the C++ interface of CRootBox. \
                   Install Python 3 and rerun CMake to build the py_rootbox library.")
endif ()
