# Let pybind11 locate and configure Python properly.
set(PYBIND11_FINDPYTHON ON)

# Allow building CPlantBox as a shared library to respect LGPL linking expectations.
option(BUILD_CPLANTBOX_SHARED "Build CPlantBox as a shared library (DLL/.so) when possible" ON)

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


if (BUILD_CPLANTBOX_SHARED)
  # Shared library (DLL on Windows, .so/.dylib on Unix)
  add_library(CPlantBox SHARED ${CPLANTBOX_SOURCES})
  if (MSVC)
    # If the project doesn't use explicit __declspec(dllexport) annotations
    # ensure the linker produces an import library by exporting all symbols.
    set_target_properties(CPlantBox PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS ON)
  endif()
else()
  # Static library containing the core C++ logic.
  add_library(CPlantBox STATIC ${CPLANTBOX_SOURCES})
endif()
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
  PUBLIC tinyxml2::tinyxml2
  PRIVATE Eigen3::Eigen
)

# Python extension that embeds CPlantBox.
pybind11_add_module(plantbox SHARED PyPlantBox.cpp)
target_compile_features(plantbox PUBLIC cxx_std_17)
target_link_libraries(plantbox PRIVATE CPlantBox)
add_dependencies(plantbox CPlantBox)

# Multi-config output locations (especially for MSVC).
if (MSVC)
  if (BUILD_CPLANTBOX_SHARED)
    set_target_properties(CPlantBox PROPERTIES
      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/$<CONFIG>"
      LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/$<CONFIG>"
      ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/$<CONFIG>"
    )
  else()
    # Static: archives (.lib) go to the top-level archive dir too
    set_target_properties(CPlantBox PROPERTIES
      ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/$<CONFIG>"
    )
  endif()
  set_target_properties(plantbox PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY_RELEASE "${CMAKE_BINARY_DIR}/Release"
    RUNTIME_OUTPUT_DIRECTORY_DEBUG "${CMAKE_BINARY_DIR}/Debug"
  )
endif()

set(CPLANTBOX_INCLUDE_DIRS
  "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/structural"
  "${CMAKE_CURRENT_SOURCE_DIR}/functional"
  "${CMAKE_CURRENT_SOURCE_DIR}/visualisation"
  "${CMAKE_CURRENT_SOURCE_DIR}/external"
)
set(CPLANTBOX_LIBRARY_DIR "${CMAKE_CURRENT_BINARY_DIR}")
set(CPLANTBOX_VERSION "2.0w")

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/CPlantBoxConfigVersion.cmake"
  VERSION ${CPLANTBOX_VERSION}
  COMPATIBILITY AnyNewerVersion
)

export(TARGETS CPlantBox
  FILE "${CMAKE_CURRENT_BINARY_DIR}/CPlantBoxTargets.cmake"
  NAMESPACE CPlantBox::
)
configure_package_config_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/../cmake/CPlantBoxConfig.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/CPlantBoxConfig.cmake"
  INSTALL_DESTINATION "lib/cmake/CPlantBox"
)
install(TARGETS CPlantBox
  EXPORT CPlantBoxTargets
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  INCLUDES DESTINATION include
)

install(EXPORT CPlantBoxTargets
  FILE CPlantBoxTargets.cmake
  NAMESPACE CPlantBox::
  DESTINATION lib/cmake/CPlantBox
)
install(FILES
  "${CMAKE_CURRENT_BINARY_DIR}/CPlantBoxConfig.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/CPlantBoxConfigVersion.cmake"
  DESTINATION lib/cmake/CPlantBox
)
install(FILES
  "${CMAKE_CURRENT_BINARY_DIR}/CPlantBoxConfig.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/CPlantBoxConfigVersion.cmake"
  DESTINATION lib/cmake/CPlantBox
)
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/"
  DESTINATION include/CPlantBox
  FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp"
)
