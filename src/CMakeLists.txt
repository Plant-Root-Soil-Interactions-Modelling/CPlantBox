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
include_directories(${PROJECT_SOURCE_DIR}/src/external/suitsparse/include/) # for PiafMunch (fallback when bundled)
include_directories(${PROJECT_SOURCE_DIR}/src/external/sundials/include/) # for PiafMunch (fallback when bundled)


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


# Dependency strategy options (Phase 4)
option(USE_SYSTEM_SUITESPARSE "Use system-installed SuiteSparse" OFF)
option(USE_SYSTEM_SUNDIALS "Use system-installed SUNDIALS" OFF)
option(BUNDLE_SUITESPARSE "Bundle vendored SuiteSparse static libraries" ON)
option(BUNDLE_SUNDIALS "Bundle vendored SUNDIALS static libraries" ON)

# Prevent conflicting selections
if(USE_SYSTEM_SUITESPARSE AND BUNDLE_SUITESPARSE)
  message(FATAL_ERROR "Both USE_SYSTEM_SUITESPARSE and BUNDLE_SUITESPARSE are ON. Choose one.")
endif()
if(USE_SYSTEM_SUNDIALS AND BUNDLE_SUNDIALS)
  message(FATAL_ERROR "Both USE_SYSTEM_SUNDIALS and BUNDLE_SUNDIALS are ON. Choose one.")
endif()

# For wheel builds (SKBUILD), default to bundled libs unless explicitly overridden
if(SKBUILD)
  if(NOT DEFINED CACHE{USE_SYSTEM_SUITESPARSE})
    set(USE_SYSTEM_SUITESPARSE OFF CACHE BOOL "" FORCE)
    set(BUNDLE_SUITESPARSE ON CACHE BOOL "" FORCE)
  endif()
  if(NOT DEFINED CACHE{USE_SYSTEM_SUNDIALS})
    set(USE_SYSTEM_SUNDIALS OFF CACHE BOOL "" FORCE)
    set(BUNDLE_SUNDIALS ON CACHE BOOL "" FORCE)
  endif()
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
add_library(libamd SHARED IMPORTED)
add_library(libbtf SHARED IMPORTED)
add_library(libcolamd SHARED IMPORTED)
add_library(suitesparseconfig SHARED IMPORTED)
if(USE_SYSTEM_SUITESPARSE)
  find_library(SS_KLU klu)
  find_library(SS_AMD amd)
  find_library(SS_COLAMD colamd)
  find_library(SS_BTF btf)
  find_library(SS_CONFIG suitesparseconfig)
  if(NOT SS_KLU OR NOT SS_AMD OR NOT SS_COLAMD OR NOT SS_BTF OR NOT SS_CONFIG)
    message(FATAL_ERROR "Requested USE_SYSTEM_SUITESPARSE but required libraries not found (klu, amd, colamd, btf, suitesparseconfig)")
  endif()
  set_property(TARGET libklu PROPERTY IMPORTED_LOCATION "${SS_KLU}")
  set_property(TARGET libamd PROPERTY IMPORTED_LOCATION "${SS_AMD}")
  set_property(TARGET libcolamd PROPERTY IMPORTED_LOCATION "${SS_COLAMD}")
  set_property(TARGET libbtf PROPERTY IMPORTED_LOCATION "${SS_BTF}")
  set_property(TARGET suitesparseconfig PROPERTY IMPORTED_LOCATION "${SS_CONFIG}")
else()
  _set_imported(libklu "${SUITESPARSE_LIB_BIN}/libklu.a" "${SUITESPARSE_LIB_SRC}/libklu.a")
  _set_imported(libamd "${SUITESPARSE_LIB_BIN}/libamd.a" "${SUITESPARSE_LIB_SRC}/libamd.a")
  _set_imported(libbtf "${SUITESPARSE_LIB_BIN}/libbtf.a" "${SUITESPARSE_LIB_SRC}/libbtf.a")
  _set_imported(libcolamd "${SUITESPARSE_LIB_BIN}/libcolamd.a" "${SUITESPARSE_LIB_SRC}/libcolamd.a")
  _set_imported(suitesparseconfig "${SUITESPARSE_LIB_BIN}/libsuitesparseconfig.a" "${SUITESPARSE_LIB_SRC}/libsuitesparseconfig.a")
endif()

add_library(arkode SHARED IMPORTED)
add_library(sundials_cvode SHARED IMPORTED)
add_library(sundials_nvecserial SHARED IMPORTED)
add_library(sundials_sunlinsolband SHARED IMPORTED)
add_library(sundials_sunlinsoldense SHARED IMPORTED)
add_library(sundials_sunlinsolklu SHARED IMPORTED)
add_library(sundials_sunlinsolpcg SHARED IMPORTED)
add_library(sundials_sunlinsolspbcgs SHARED IMPORTED)
add_library(sundials_sunlinsolspfgmr SHARED IMPORTED)
add_library(sundials_sunlinsolspgmr SHARED IMPORTED)
add_library(sundials_sunlinsolsptfqmr SHARED IMPORTED)
add_library(sundials_sunmatrixband SHARED IMPORTED)
add_library(sundials_sunmatrixdense SHARED IMPORTED)
add_library(sundials_sunmatrixsparse SHARED IMPORTED)
add_library(sundials_sunnonlinsolnewton SHARED IMPORTED)
if(USE_SYSTEM_SUNDIALS)
  find_library(SD_ARKODE sundials_arkode)
  find_library(SD_CVODE sundials_cvode)
  find_library(SD_NVECSERIAL sundials_nvecserial)
  find_library(SD_SUNLINSOLBAND sundials_sunlinsolband)
  find_library(SD_SUNLINSOLDENSE sundials_sunlinsoldense)
  find_library(SD_SUNLINSOLKLU sundials_sunlinsolklu)
  find_library(SD_SUNLINSOLPCG sundials_sunlinsolpcg)
  find_library(SD_SUNLINSOLSPBCGS sundials_sunlinsolspbcgs)
  find_library(SD_SUNLINSOLSPFGMR sundials_sunlinsolspfgmr)
  find_library(SD_SUNLINSOLSPGMR sundials_sunlinsolspgmr)
  find_library(SD_SUNLINSOLSPTFQMR sundials_sunlinsolsptfqmr)
  find_library(SD_SUNMATRIXBAND sundials_sunmatrixband)
  find_library(SD_SUNMATRIXDENSE sundials_sunmatrixdense)
  find_library(SD_SUNMATRIXSPARSE sundials_sunmatrixsparse)
  find_library(SD_SUNNONLINSOLNEWTON sundials_sunnonlinsolnewton)
  if(NOT SD_CVODE OR NOT SD_NVECSERIAL)
    message(FATAL_ERROR "Requested USE_SYSTEM_SUNDIALS but required libraries not fully found")
  endif()
  if(SD_ARKODE)
    set_property(TARGET arkode PROPERTY IMPORTED_LOCATION "${SD_ARKODE}")
  endif()
  set_property(TARGET sundials_cvode PROPERTY IMPORTED_LOCATION "${SD_CVODE}")
  set_property(TARGET sundials_nvecserial PROPERTY IMPORTED_LOCATION "${SD_NVECSERIAL}")
  if(SD_SUNLINSOLBAND)
    set_property(TARGET sundials_sunlinsolband PROPERTY IMPORTED_LOCATION "${SD_SUNLINSOLBAND}")
  endif()
  if(SD_SUNLINSOLDENSE)
    set_property(TARGET sundials_sunlinsoldense PROPERTY IMPORTED_LOCATION "${SD_SUNLINSOLDENSE}")
  endif()
  if(SD_SUNLINSOLKLU)
    set_property(TARGET sundials_sunlinsolklu PROPERTY IMPORTED_LOCATION "${SD_SUNLINSOLKLU}")
  endif()
  if(SD_SUNLINSOLPCG)
    set_property(TARGET sundials_sunlinsolpcg PROPERTY IMPORTED_LOCATION "${SD_SUNLINSOLPCG}")
  endif()
  if(SD_SUNLINSOLSPBCGS)
    set_property(TARGET sundials_sunlinsolspbcgs PROPERTY IMPORTED_LOCATION "${SD_SUNLINSOLSPBCGS}")
  endif()
  if(SD_SUNLINSOLSPFGMR)
    set_property(TARGET sundials_sunlinsolspfgmr PROPERTY IMPORTED_LOCATION "${SD_SUNLINSOLSPFGMR}")
  endif()
  if(SD_SUNLINSOLSPGMR)
    set_property(TARGET sundials_sunlinsolspgmr PROPERTY IMPORTED_LOCATION "${SD_SUNLINSOLSPGMR}")
  endif()
  if(SD_SUNLINSOLSPTFQMR)
    set_property(TARGET sundials_sunlinsolsptfqmr PROPERTY IMPORTED_LOCATION "${SD_SUNLINSOLSPTFQMR}")
  endif()
  if(SD_SUNMATRIXBAND)
    set_property(TARGET sundials_sunmatrixband PROPERTY IMPORTED_LOCATION "${SD_SUNMATRIXBAND}")
  endif()
  if(SD_SUNMATRIXDENSE)
    set_property(TARGET sundials_sunmatrixdense PROPERTY IMPORTED_LOCATION "${SD_SUNMATRIXDENSE}")
  endif()
  if(SD_SUNMATRIXSPARSE)
    set_property(TARGET sundials_sunmatrixsparse PROPERTY IMPORTED_LOCATION "${SD_SUNMATRIXSPARSE}")
  endif()
  if(SD_SUNNONLINSOLNEWTON)
    set_property(TARGET sundials_sunnonlinsolnewton PROPERTY IMPORTED_LOCATION "${SD_SUNNONLINSOLNEWTON}")
  endif()
else()
  _set_imported(arkode "${SUNDIALS_LIB_BIN}/libsundials_arkode.a" "${SUNDIALS_LIB_SRC}/libsundials_arkode.a")
  _set_imported(sundials_cvode "${SUNDIALS_LIB_BIN}/libsundials_cvode.a" "${SUNDIALS_LIB_SRC}/libsundials_cvode.a")
  _set_imported(sundials_nvecserial "${SUNDIALS_LIB_BIN}/libsundials_nvecserial.a" "${SUNDIALS_LIB_SRC}/libsundials_nvecserial.a")
  _set_imported(sundials_sunlinsolband "${SUNDIALS_LIB_BIN}/libsundials_sunlinsolband.a" "${SUNDIALS_LIB_SRC}/libsundials_sunlinsolband.a")
  _set_imported(sundials_sunlinsoldense "${SUNDIALS_LIB_BIN}/libsundials_sunlinsoldense.a" "${SUNDIALS_LIB_SRC}/libsundials_sunlinsoldense.a")
  _set_imported(sundials_sunlinsolklu "${SUNDIALS_LIB_BIN}/libsundials_sunlinsolklu.a" "${SUNDIALS_LIB_SRC}/libsundials_sunlinsolklu.a")
  _set_imported(sundials_sunlinsolpcg "${SUNDIALS_LIB_BIN}/libsundials_sunlinsolpcg.a" "${SUNDIALS_LIB_SRC}/libsundials_sunlinsolpcg.a")
  _set_imported(sundials_sunlinsolspbcgs "${SUNDIALS_LIB_BIN}/libsundials_sunlinsolspbcgs.a" "${SUNDIALS_LIB_SRC}/libsundials_sunlinsolspbcgs.a")
  _set_imported(sundials_sunlinsolspfgmr "${SUNDIALS_LIB_BIN}/libsundials_sunlinsolspfgmr.a" "${SUNDIALS_LIB_SRC}/libsundials_sunlinsolspfgmr.a")
  _set_imported(sundials_sunlinsolspgmr "${SUNDIALS_LIB_BIN}/libsundials_sunlinsolspgmr.a" "${SUNDIALS_LIB_SRC}/libsundials_sunlinsolspgmr.a")
  _set_imported(sundials_sunlinsolsptfqmr "${SUNDIALS_LIB_BIN}/libsundials_sunlinsolsptfqmr.a" "${SUNDIALS_LIB_SRC}/libsundials_sunlinsolsptfqmr.a")
  _set_imported(sundials_sunmatrixband "${SUNDIALS_LIB_BIN}/libsundials_sunmatrixband.a" "${SUNDIALS_LIB_SRC}/libsundials_sunmatrixband.a")
  _set_imported(sundials_sunmatrixdense "${SUNDIALS_LIB_BIN}/libsundials_sunmatrixdense.a" "${SUNDIALS_LIB_SRC}/libsundials_sunmatrixdense.a")
  _set_imported(sundials_sunmatrixsparse "${SUNDIALS_LIB_BIN}/libsundials_sunmatrixsparse.a" "${SUNDIALS_LIB_SRC}/libsundials_sunmatrixsparse.a")
  _set_imported(sundials_sunnonlinsolnewton "${SUNDIALS_LIB_BIN}/libsundials_sunnonlinsolnewton.a" "${SUNDIALS_LIB_SRC}/libsundials_sunnonlinsolnewton.a")
endif()
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
set(SS_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/src/external/suitsparse/include")
set(SD_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/src/external/sundials/include")
if(USE_SYSTEM_SUITESPARSE)
  find_path(SS_INCLUDE_DIR NAMES klu.h SuiteSparse_config.h)
endif()
if(USE_SYSTEM_SUNDIALS)
  find_path(SD_INCLUDE_DIR NAMES sundials/sundials_types.h sundials/sundials_config.h)
endif()

target_include_directories(CPlantBox PUBLIC ${SS_INCLUDE_DIR} ${SD_INCLUDE_DIR})
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
# Workaround: some pybind11 tooling expects Python3::Python even for module-only builds.
# If Embed component is unavailable, provide a dummy INTERFACE target so python3_add_library does not fail.
if(NOT TARGET Python3::Python)
  add_library(Python3::Python INTERFACE IMPORTED)
endif()
if(NOT EXISTS "${CMAKE_SOURCE_DIR}/src/external/pybind11/CMakeLists.txt")
  message(FATAL_ERROR "Vendored pybind11 missing. Run: git submodule update --init --recursive")
endif()
add_subdirectory(external/pybind11)

pybind11_add_module(_plantbox SHARED 
            PyPlantBox.cpp                                  
            )
# Link Python module; on Linux wheel builds, embed static CPlantBox objects with whole-archive
if(SKBUILD AND UNIX AND NOT APPLE)
  target_link_libraries(_plantbox PRIVATE Python3::Module "-Wl,--whole-archive" CPlantBox "-Wl,--no-whole-archive")
else()
  target_link_libraries(_plantbox PUBLIC CPlantBox Python3::Module)
endif()

# Ensure static CPlantBox is fully included in the extension on Linux wheel builds


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






