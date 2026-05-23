# -----------------------------
# Build options
# -----------------------------
option(CPB_BUILD_PYTHON "Build the Python bindings" ON)
option(CPB_BUILD_SHARED "Build CPlantBox as a shared library" ON)

set(CPB_SUITESPARSE_PROVIDER "bundled" CACHE STRING "SuiteSparse provider: bundled, source, or system")
set_property(CACHE CPB_SUITESPARSE_PROVIDER PROPERTY STRINGS bundled source system)
set(CPB_SUNDIALS_PROVIDER "bundled" CACHE STRING "SUNDIALS provider: bundled, source, or system")
set_property(CACHE CPB_SUNDIALS_PROVIDER PROPERTY STRINGS bundled source system)
set(CPB_NATIVE_DEPS_PREFIX "" CACHE PATH "Prefix containing source-built SuiteSparse/SUNDIALS dependencies")

string(TOLOWER "${CPB_SUITESPARSE_PROVIDER}" CPB_SUITESPARSE_PROVIDER)
string(TOLOWER "${CPB_SUNDIALS_PROVIDER}" CPB_SUNDIALS_PROVIDER)
set(CPB_SUITESPARSE_PROVIDER "${CPB_SUITESPARSE_PROVIDER}" CACHE STRING "SuiteSparse provider: bundled, source, or system" FORCE)
set(CPB_SUNDIALS_PROVIDER "${CPB_SUNDIALS_PROVIDER}" CACHE STRING "SUNDIALS provider: bundled, source, or system" FORCE)

message(STATUS "Project dir: ${PROJECT_SOURCE_DIR}")
message(STATUS "Current source dir: ${CMAKE_CURRENT_SOURCE_DIR}")
message(STATUS "Project binary dir: ${PROJECT_BINARY_DIR}")

if(APPLE)
    set(CPB_INSTALL_RPATH "@loader_path")
else()
    set(CPB_INSTALL_RPATH "$ORIGIN")
endif()

# -----------------------------
# Core library
# -----------------------------
if(CPB_BUILD_SHARED)
    set(CPB_LIBRARY_TYPE SHARED)
else()
    set(CPB_LIBRARY_TYPE STATIC)
endif()

add_library(CPlantBox ${CPB_LIBRARY_TYPE})
target_sources(CPlantBox
    PRIVATE
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
)

set_target_properties(CPlantBox PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
    BUILD_WITH_INSTALL_RPATH TRUE
    INSTALL_RPATH "${CPB_INSTALL_RPATH}"
)

target_include_directories(CPlantBox
    PUBLIC
        "${CMAKE_CURRENT_SOURCE_DIR}/structural"
        "${CMAKE_CURRENT_SOURCE_DIR}/functional"
        "${CMAKE_CURRENT_SOURCE_DIR}/visualisation"
        "${CMAKE_CURRENT_SOURCE_DIR}/external/tinyxml2"
        "${CMAKE_CURRENT_SOURCE_DIR}/external/eigen"
        "${CMAKE_CURRENT_SOURCE_DIR}/external"
)

# -----------------------------
# Optional features
# -----------------------------
if(PIAFMUNCH)
    include("${CMAKE_CURRENT_LIST_DIR}/../cmake/CPlantBoxPiafMunch.cmake")
    cplantbox_configure_piafmunch(CPlantBox)
endif()

if(CPB_BUILD_PYTHON)
    include("${CMAKE_CURRENT_LIST_DIR}/../cmake/CPlantBoxPython.cmake")
    cplantbox_configure_python_bindings(CPlantBox)
else()
    include(GNUInstallDirs)
    install(TARGETS CPlantBox
            LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
            ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
            RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}")
endif()
