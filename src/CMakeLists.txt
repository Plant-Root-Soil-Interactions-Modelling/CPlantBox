# Let pybind11 locate and configure Python properly.
set(PYBIND11_FINDPYTHON ON)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(tinyxml2 CONFIG REQUIRED)
find_package(pybind11 CONFIG REQUIRED)

add_compile_definitions(NOMINMAX) # otherwise causes errors with perirhizal.cpp..

set(CPLANTBOX_SOURCES
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
  external/aabbcc/AABB.cc
  external/gauss_legendre/gauss_legendre.cpp
)
# Should be all except for PiafMunch, legacy SUNDIALS makes problems with vcpkg...


# Static library containing the core C++ logic.
add_library(CPlantBox STATIC ${CPLANTBOX_SOURCES})
target_compile_features(CPlantBox PUBLIC cxx_std_17)

target_include_directories(CPlantBox
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/structural>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/functional>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/visualisation>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/external>
)

target_link_libraries(CPlantBox
  PRIVATE tinyxml2::tinyxml2
  PRIVATE Eigen3::Eigen
)

# Python extension that embeds CPlantBox.
pybind11_add_module(plantbox SHARED PyPlantBox.cpp)
target_compile_features(plantbox PUBLIC cxx_std_17)
target_link_libraries(plantbox PRIVATE CPlantBox)
add_dependencies(plantbox CPlantBox)

# Multi-config output locations (especially for MSVC).
if (MSVC)
  set_target_properties(CPlantBox PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY_RELEASE "${CMAKE_BINARY_DIR}/Release"
    ARCHIVE_OUTPUT_DIRECTORY_DEBUG "${CMAKE_BINARY_DIR}/Debug"
  )
  set_target_properties(plantbox PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY_RELEASE "${CMAKE_BINARY_DIR}/Release"
    RUNTIME_OUTPUT_DIRECTORY_DEBUG "${CMAKE_BINARY_DIR}/Debug"
  )
endif()
