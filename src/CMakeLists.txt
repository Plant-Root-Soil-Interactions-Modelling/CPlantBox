cmake_minimum_required(VERSION 3.16)
project(CPlantBox LANGUAGES CXX)

message(STATUS "Current source dir: ${CMAKE_CURRENT_SOURCE_DIR}")
message(STATUS "Project binary dir: ${PROJECT_BINARY_DIR}")

# 1. Dependencies
find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
add_subdirectory(external/pybind11)

# 2. Main C++ library
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

add_library(libklu SHARED IMPORTED)
set_property(TARGET libklu PROPERTY IMPORTED_LOCATION ${PROJECT_BINARY_DIR}/external/suitsparse/lib/libklu.a)
add_library(libamd SHARED IMPORTED)
set_property(TARGET libamd PROPERTY IMPORTED_LOCATION ${PROJECT_BINARY_DIR}/external/suitsparse/lib/libamd.a)
add_library(libbtf SHARED IMPORTED)
set_property(TARGET libbtf PROPERTY IMPORTED_LOCATION ${PROJECT_BINARY_DIR}/external/suitsparse/lib/libbtf.a)
add_library(libcolamd SHARED IMPORTED)
set_property(TARGET libcolamd PROPERTY IMPORTED_LOCATION ${PROJECT_BINARY_DIR}/external/suitsparse/lib/libcolamd.a)
add_library(suitesparseconfig SHARED IMPORTED)
set_property(TARGET suitesparseconfig PROPERTY IMPORTED_LOCATION ${PROJECT_BINARY_DIR}/external/suitsparse/lib/libsuitesparseconfig.a)
add_library(arkode SHARED IMPORTED)
set_target_properties(arkode PROPERTIES IMPORTED_LOCATION ${PROJECT_BINARY_DIR}/external/sundials/lib/libsundials_arkode.a)
add_library(sundials_cvode SHARED IMPORTED)
set_property(TARGET sundials_cvode PROPERTY IMPORTED_LOCATION ${PROJECT_BINARY_DIR}/external/sundials/lib/libsundials_cvode.a)
add_library(sundials_nvecserial SHARED IMPORTED)
set_property(TARGET sundials_nvecserial PROPERTY IMPORTED_LOCATION ${PROJECT_BINARY_DIR}/external/sundials/lib/libsundials_nvecserial.a)
add_library(sundials_sunlinsolband SHARED IMPORTED)
set_property(TARGET sundials_sunlinsolband PROPERTY IMPORTED_LOCATION ${PROJECT_BINARY_DIR}/external/sundials/lib/libsundials_sunlinsolband.a)
add_library(sundials_sunlinsoldense SHARED IMPORTED)
set_property(TARGET sundials_sunlinsoldense PROPERTY IMPORTED_LOCATION ${PROJECT_BINARY_DIR}/external/sundials/lib/libsundials_sunlinsoldense.a)
add_library(sundials_sunlinsolklu SHARED IMPORTED)
set_property(TARGET sundials_sunlinsolklu PROPERTY IMPORTED_LOCATION ${PROJECT_BINARY_DIR}/external/sundials/lib/libsundials_sunlinsolklu.a)
add_library(sundials_sunlinsolpcg SHARED IMPORTED)
set_property(TARGET sundials_sunlinsolpcg PROPERTY IMPORTED_LOCATION ${PROJECT_BINARY_DIR}/external/sundials/lib/libsundials_sunlinsolpcg.a)
add_library(sundials_sunlinsolspbcgs SHARED IMPORTED)
set_property(TARGET sundials_sunlinsolspbcgs PROPERTY IMPORTED_LOCATION ${PROJECT_BINARY_DIR}/external/sundials/lib/libsundials_sunlinsolspbcgs.a)
add_library(sundials_sunlinsolspfgmr SHARED IMPORTED)
set_property(TARGET sundials_sunlinsolspfgmr PROPERTY IMPORTED_LOCATION ${PROJECT_BINARY_DIR}/external/sundials/lib/libsundials_sunlinsolspfgmr.a)
add_library(sundials_sunlinsolspgmr SHARED IMPORTED)
set_property(TARGET sundials_sunlinsolspgmr PROPERTY IMPORTED_LOCATION ${PROJECT_BINARY_DIR}/external/sundials/lib/libsundials_sunlinsolspgmr.a)
add_library(sundials_sunlinsolsptfqmr SHARED IMPORTED)
set_property(TARGET sundials_sunlinsolsptfqmr PROPERTY IMPORTED_LOCATION ${PROJECT_BINARY_DIR}/external/sundials/lib/libsundials_sunlinsolsptfqmr.a)
add_library(sundials_sunmatrixband SHARED IMPORTED)
set_property(TARGET sundials_sunmatrixband PROPERTY IMPORTED_LOCATION ${PROJECT_BINARY_DIR}/external/sundials/lib/libsundials_sunmatrixband.a)
add_library(sundials_sunmatrixdense SHARED IMPORTED)
set_property(TARGET sundials_sunmatrixdense PROPERTY IMPORTED_LOCATION ${PROJECT_BINARY_DIR}/external/sundials/lib/libsundials_sunmatrixdense.a)
add_library(sundials_sunmatrixsparse SHARED IMPORTED)
set_target_properties(sundials_sunmatrixsparse PROPERTIES IMPORTED_LOCATION ${PROJECT_BINARY_DIR}/external/sundials/lib/libsundials_sunmatrixsparse.a)
add_library(sundials_sunnonlinsolnewton SHARED IMPORTED)
set_target_properties(sundials_sunnonlinsolnewton PROPERTIES IMPORTED_LOCATION ${PROJECT_BINARY_DIR}/external/sundials/lib/libsundials_sunnonlinsolnewton.a)
target_link_libraries(CPlantBox 
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
						suitesparseconfig
						arkode 
						libbtf
						libcolamd 
						libamd 
						libklu
						)

set_target_properties(CPlantBox PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/)
target_include_directories(
                        CPlantBox
                        PUBLIC
                        ${PROJECT_BINARY_DIR}/external/sundials/include
                        ${PROJECT_BINARY_DIR}/external/suitsparse/include
)

# Set library output directory
set_target_properties(CPlantBox PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
    BUILD_WITH_INSTALL_RPATH TRUE
    INSTALL_RPATH "$ORIGIN"
)


# -----------------------------
# 6. Python binding
# -----------------------------

pybind11_add_module(plantbox PyPlantBox.cpp)
target_link_libraries(plantbox PRIVATE CPlantBox Python3::Module)

set_target_properties(plantbox PROPERTIES
    INSTALL_RPATH "$ORIGIN"
)

# Install plantbox module into Python site-packages
install(TARGETS CPlantBox LIBRARY DESTINATION "${Python3_SITEARCH}")
install(TARGETS plantbox LIBRARY DESTINATION "${Python3_SITEARCH}")
