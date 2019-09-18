#
# 1. Make RootBox library
#

# file(GLOB all_sources *.cpp) # i changed that to exclude PythonRootSystem.cpp
add_library(CRootBox
            organparameter.cpp            
            rootparameter.cpp
			seedparameter.cpp  
			leafparameter.cpp          
            Organ.cpp            
            Root.cpp
            Seed.cpp
            Organism.cpp
            RootSystem.cpp
            analysis.cpp
            sdf.cpp
            tropism.cpp
			../external/tinyxml2.cpp            
            ../external/aabbcc/AABB.cc)
            #../external/gauss_legendre/gauss_legendre.cpp
			

set_target_properties(CRootBox PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)

#
# 2. Make py_rootbox library
#

# Find necessary packages for py_rootbox
find_package(PythonLibs 3)
find_package(PythonInterp 3)    # creates PYTHON_VERSION_MAJOR, PYTHON_VERSION_MINOR

# Only build the library if Python 3 was found
if (PythonLibs_FOUND)
  include_directories( ${PYTHON_INCLUDE_DIRS} )
#
#  find_package( Boost COMPONENTS python${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR} REQUIRED )
	find_package( Boost COMPONENTS python-py${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR} REQUIRED ) # if above is not working
  include_directories( ${Boost_INCLUDE_DIR} )

  # make shared library for Pyhton binding
  add_library(py_rootbox SHARED
            organparameter.cpp            
            rootparameter.cpp
			seedparameter.cpp      
			leafparameter.cpp      
			seedparameter.cpp
			stemparameter.cpp
			#plant.cpp
			Stem.cpp
            Organ.cpp            
            Root.cpp
            Seed.cpp
            Organism.cpp
            RootSystem.cpp
            PythonRootSystem.cpp            
            analysis.cpp
            sdf.cpp
            tropism.cpp
			../external/tinyxml2/tinyxml2.cpp                 
            ../external/aabbcc/AABB.cc
            #../external/gauss_legendre/gauss_legendre.cpp
)
  target_link_libraries(py_rootbox ${Boost_LIBRARIES} )
  set_target_properties(py_rootbox PROPERTIES PREFIX "" )
  set_target_properties(py_rootbox PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/python)
else ()
  message(WARNING "Python 3 was not found on your system! You can only use the C++ interface of CRootBox. \
                   Install Python 3 and rerun CMake to build the py_rootbox library.")
endif ()
