# Disable LTO/IPO to avoid toolchain issues in container builds
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION OFF)
#
# 1. Make CPlantBox library
#

include_directories(structural)
include_directories(functional)
include_directories(visualisation)

include_directories(external/tinyxml2)
include_directories(external/eigen) # fetching as git submodule @tag 3.3.7

include_directories(external/) # for PiafMunch
include_directories(${PROJECT_SOURCE_DIR}/src/external/suitsparse/include/) # for PiafMunch
include_directories(${PROJECT_SOURCE_DIR}/src/external/sundials/include/) # for PiafMunch


# Build CPlantBox as STATIC when building wheels (scikit-build), to embed into the
# Python extension and avoid shipping a separate shared library in wheels.
set(CPB_BUILD_SHARED ON)
if(SKBUILD)
  set(CPB_BUILD_SHARED OFF)
endif()

if(CPB_BUILD_SHARED)
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
			functional/Photosynthesis.cpp
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
			external/PiafMunch/runPM.cpp     
			external/PiafMunch/PiafMunch2.cpp     
			external/PiafMunch/index_vector.cpp
			external/PiafMunch/initialize.cpp
			external/PiafMunch/odepack.cpp
			external/PiafMunch/PM_KLU.cpp
			external/PiafMunch/PM_matrix.cpp
			external/PiafMunch/PM_vector.cpp
			external/PiafMunch/solve.cpp
            external/PiafMunch/sparse_matrix.cpp            
)
else()
add_library(CPlantBox STATIC
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
            functional/Photosynthesis.cpp
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
            external/PiafMunch/runPM.cpp     
            external/PiafMunch/PiafMunch2.cpp     
            external/PiafMunch/index_vector.cpp
            external/PiafMunch/initialize.cpp
            external/PiafMunch/odepack.cpp
            external/PiafMunch/PM_KLU.cpp
            external/PiafMunch/PM_matrix.cpp
            external/PiafMunch/PM_vector.cpp
            external/PiafMunch/solve.cpp
            external/PiafMunch/sparse_matrix.cpp            
)
set_target_properties(CPlantBox PROPERTIES POSITION_INDEPENDENT_CODE ON)
endif()


set(SUITESPARSE_LIB_BIN "${PROJECT_BINARY_DIR}/src/external/suitsparse/lib")
set(SUITESPARSE_LIB_SRC "${CMAKE_SOURCE_DIR}/src/external/suitsparse/lib")
set(SUNDIALS_LIB_BIN   "${PROJECT_BINARY_DIR}/src/external/sundials/lib")
set(SUNDIALS_LIB_SRC   "${CMAKE_SOURCE_DIR}/src/external/sundials/lib")

function(_set_imported tgt bin_path src_path)
    if(EXISTS "${bin_path}")
        set_property(TARGET ${tgt} PROPERTY IMPORTED_LOCATION "${bin_path}")
    elseif(EXISTS "${src_path}")
        set_property(TARGET ${tgt} PROPERTY IMPORTED_LOCATION "${src_path}")
    else()
        message(FATAL_ERROR "Missing archive for ${tgt}: ${bin_path} and ${src_path} not found")
    endif()
endfunction()

add_library(libklu SHARED IMPORTED)
_set_imported(libklu "${SUITESPARSE_LIB_BIN}/libklu.a" "${SUITESPARSE_LIB_SRC}/libklu.a")
add_library(libamd SHARED IMPORTED)
_set_imported(libamd "${SUITESPARSE_LIB_BIN}/libamd.a" "${SUITESPARSE_LIB_SRC}/libamd.a")
add_library(libbtf SHARED IMPORTED)
_set_imported(libbtf "${SUITESPARSE_LIB_BIN}/libbtf.a" "${SUITESPARSE_LIB_SRC}/libbtf.a")
add_library(libcolamd SHARED IMPORTED)
_set_imported(libcolamd "${SUITESPARSE_LIB_BIN}/libcolamd.a" "${SUITESPARSE_LIB_SRC}/libcolamd.a")
add_library(suitesparseconfig SHARED IMPORTED)
_set_imported(suitesparseconfig "${SUITESPARSE_LIB_BIN}/libsuitesparseconfig.a" "${SUITESPARSE_LIB_SRC}/libsuitesparseconfig.a")

add_library(arkode SHARED IMPORTED)
_set_imported(arkode "${SUNDIALS_LIB_BIN}/libsundials_arkode.a" "${SUNDIALS_LIB_SRC}/libsundials_arkode.a")
add_library(sundials_cvode SHARED IMPORTED)
_set_imported(sundials_cvode "${SUNDIALS_LIB_BIN}/libsundials_cvode.a" "${SUNDIALS_LIB_SRC}/libsundials_cvode.a")
add_library(sundials_nvecserial SHARED IMPORTED)
_set_imported(sundials_nvecserial "${SUNDIALS_LIB_BIN}/libsundials_nvecserial.a" "${SUNDIALS_LIB_SRC}/libsundials_nvecserial.a")
add_library(sundials_sunlinsolband SHARED IMPORTED)
_set_imported(sundials_sunlinsolband "${SUNDIALS_LIB_BIN}/libsundials_sunlinsolband.a" "${SUNDIALS_LIB_SRC}/libsundials_sunlinsolband.a")
add_library(sundials_sunlinsoldense SHARED IMPORTED)
_set_imported(sundials_sunlinsoldense "${SUNDIALS_LIB_BIN}/libsundials_sunlinsoldense.a" "${SUNDIALS_LIB_SRC}/libsundials_sunlinsoldense.a")
add_library(sundials_sunlinsolklu SHARED IMPORTED)
_set_imported(sundials_sunlinsolklu "${SUNDIALS_LIB_BIN}/libsundials_sunlinsolklu.a" "${SUNDIALS_LIB_SRC}/libsundials_sunlinsolklu.a")
add_library(sundials_sunlinsolpcg SHARED IMPORTED)
_set_imported(sundials_sunlinsolpcg "${SUNDIALS_LIB_BIN}/libsundials_sunlinsolpcg.a" "${SUNDIALS_LIB_SRC}/libsundials_sunlinsolpcg.a")
add_library(sundials_sunlinsolspbcgs SHARED IMPORTED)
_set_imported(sundials_sunlinsolspbcgs "${SUNDIALS_LIB_BIN}/libsundials_sunlinsolspbcgs.a" "${SUNDIALS_LIB_SRC}/libsundials_sunlinsolspbcgs.a")
add_library(sundials_sunlinsolspfgmr SHARED IMPORTED)
_set_imported(sundials_sunlinsolspfgmr "${SUNDIALS_LIB_BIN}/libsundials_sunlinsolspfgmr.a" "${SUNDIALS_LIB_SRC}/libsundials_sunlinsolspfgmr.a")
add_library(sundials_sunlinsolspgmr SHARED IMPORTED)
_set_imported(sundials_sunlinsolspgmr "${SUNDIALS_LIB_BIN}/libsundials_sunlinsolspgmr.a" "${SUNDIALS_LIB_SRC}/libsundials_sunlinsolspgmr.a")
add_library(sundials_sunlinsolsptfqmr SHARED IMPORTED)
_set_imported(sundials_sunlinsolsptfqmr "${SUNDIALS_LIB_BIN}/libsundials_sunlinsolsptfqmr.a" "${SUNDIALS_LIB_SRC}/libsundials_sunlinsolsptfqmr.a")
add_library(sundials_sunmatrixband SHARED IMPORTED)
_set_imported(sundials_sunmatrixband "${SUNDIALS_LIB_BIN}/libsundials_sunmatrixband.a" "${SUNDIALS_LIB_SRC}/libsundials_sunmatrixband.a")
add_library(sundials_sunmatrixdense SHARED IMPORTED)
_set_imported(sundials_sunmatrixdense "${SUNDIALS_LIB_BIN}/libsundials_sunmatrixdense.a" "${SUNDIALS_LIB_SRC}/libsundials_sunmatrixdense.a")
add_library(sundials_sunmatrixsparse SHARED IMPORTED)
_set_imported(sundials_sunmatrixsparse "${SUNDIALS_LIB_BIN}/libsundials_sunmatrixsparse.a" "${SUNDIALS_LIB_SRC}/libsundials_sunmatrixsparse.a")
add_library(sundials_sunnonlinsolnewton SHARED IMPORTED)
_set_imported(sundials_sunnonlinsolnewton "${SUNDIALS_LIB_BIN}/libsundials_sunnonlinsolnewton.a" "${SUNDIALS_LIB_SRC}/libsundials_sunnonlinsolnewton.a")
target_link_libraries(CPlantBox PUBLIC
                        # SUNDIALS components (order mostly independent)
                        sundials_sunnonlinsolnewton 
                        sundials_sunmatrixsparse 
                        sundials_sunlinsolklu
                        sundials_sunmatrixdense 
                        sundials_sunlinsoldense
                        sundials_sunmatrixband 
                        sundials_sunlinsolband
                        sundials_sunlinsolsptfqmr 
                        sundials_nvecserial
                        sundials_sunlinsolspgmr 
                        sundials_cvode
                        sundials_sunlinsolspfgmr 
                        sundials_sunlinsolspbcgs 
                        sundials_sunlinsolpcg 
                        arkode 
                        
                        # SuiteSparse (order matters for static linking):
                        # user first, then its dependencies
                        libklu
                        libamd
                        libcolamd
                        libbtf
                        suitesparseconfig
                        )

# Avoid forcing output directories for portability; use CMake defaults (skbuild/packaging will install)
# set_target_properties(CPlantBox PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/)
target_include_directories(
                        CPlantBox
                        PUBLIC
                        ${PROJECT_BINARY_DIR}/src/external/sundials/include
                        ${PROJECT_BINARY_DIR}/src/external/suitsparse/include
)
#
# 2. Make CPlantBox Pyhthon binding
#

# Python and pybind11 (modern discovery)
# Prefer Development.Module for wheels; only require Embed outside skbuild builds
set(_CPB_PY_COMPONENTS Interpreter Development.Module)
if(NOT SKBUILD)
  list(APPEND _CPB_PY_COMPONENTS Development.Embed)
endif()
find_package(Python3 COMPONENTS ${_CPB_PY_COMPONENTS} REQUIRED)
if(NOT EXISTS "${CMAKE_SOURCE_DIR}/src/external/pybind11/CMakeLists.txt")
  message(FATAL_ERROR "Vendored pybind11 missing. Run: git submodule update --init --recursive")
endif()
add_subdirectory(external/pybind11)

pybind11_add_module(_plantbox SHARED 
            PyPlantBox.cpp                                  
            )
target_link_libraries(_plantbox PUBLIC 
            CPlantBox Python3::Module)

# Install rules to allow packaging via scikit-build-core
install(TARGETS _plantbox
        LIBRARY DESTINATION plantbox
        ARCHIVE DESTINATION plantbox
        RUNTIME DESTINATION plantbox
        )

# Install model parameter data into the wheel under plantbox/, when building via skbuild
if(SKBUILD)
  # Start conservatively: include the structural parameters used by smoke tests
  install(DIRECTORY ${CMAKE_SOURCE_DIR}/modelparameter/structural/
          DESTINATION plantbox/modelparameter/structural
          FILES_MATCHING PATTERN "*")
endif()
			
# Avoid forcing output directories; keep default so wheels/deploys are simpler
# set_target_properties(plantbox PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/)






