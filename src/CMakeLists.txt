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
# 1. Sources
# -----------------------------
set(CPB_SOURCES
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

if(PIAFMUNCH)
    list(APPEND CPB_SOURCES
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
endif()

# -----------------------------
# 2. Core library
# -----------------------------
if(CPB_BUILD_SHARED)
    set(CPB_LIBRARY_TYPE SHARED)
else()
    set(CPB_LIBRARY_TYPE STATIC)
endif()

add_library(CPlantBox ${CPB_LIBRARY_TYPE} ${CPB_SOURCES})
set_target_properties(CPlantBox PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
    BUILD_WITH_INSTALL_RPATH TRUE
    INSTALL_RPATH "${CPB_INSTALL_RPATH}"
)

target_include_directories(CPlantBox
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/structural
        ${CMAKE_CURRENT_SOURCE_DIR}/functional
        ${CMAKE_CURRENT_SOURCE_DIR}/visualisation
        ${CMAKE_CURRENT_SOURCE_DIR}/external/tinyxml2
        ${CMAKE_CURRENT_SOURCE_DIR}/external/eigen
        ${CMAKE_CURRENT_SOURCE_DIR}/external
)

# -----------------------------
# 3. Optional PiafMunch native dependencies
# -----------------------------
include("${PROJECT_SOURCE_DIR}/cmake/CPlantBoxPiafMunch.cmake")

if(PIAFMUNCH)
    cplantbox_configure_piafmunch(CPlantBox)
endif()

# -----------------------------
# 4. Python binding and package install
# -----------------------------
if(CPB_BUILD_PYTHON)
    find_package(Python3 COMPONENTS Interpreter Development.Module REQUIRED)
    message(STATUS "Python library path (site-packages): ${Python3_SITEARCH}")

    add_subdirectory(external/pybind11)

    pybind11_add_module(plantbox PyPlantBox.cpp)
    target_link_libraries(plantbox PRIVATE CPlantBox Python3::Module)

    set_target_properties(plantbox PROPERTIES
        INSTALL_RPATH "${CPB_INSTALL_RPATH}"
    )

    # Install plantbox module into Python site-packages for native CMake installs.
    # scikit-build-core stages wheel contents under the CMake install prefix, so wheel
    # builds must use a relative install destination.
    if(SKBUILD)
        set(CPLANTBOX_PYTHON_INSTALL_DIR "plantbox")
    else()
        set(CPLANTBOX_PYTHON_INSTALL_DIR "${Python3_SITEARCH}/plantbox")
    endif()

    if(CPB_BUILD_SHARED)
        install(TARGETS CPlantBox
                LIBRARY DESTINATION "${CPLANTBOX_PYTHON_INSTALL_DIR}"
                RUNTIME DESTINATION "${CPLANTBOX_PYTHON_INSTALL_DIR}")
    endif()
    install(TARGETS plantbox
            LIBRARY DESTINATION "${CPLANTBOX_PYTHON_INSTALL_DIR}"
            RUNTIME DESTINATION "${CPLANTBOX_PYTHON_INSTALL_DIR}")

    install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py DESTINATION ${CPLANTBOX_PYTHON_INSTALL_DIR})

    install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/visualisation DESTINATION ${CPLANTBOX_PYTHON_INSTALL_DIR}
            FILES_MATCHING
            PATTERN "*.py"
            PATTERN "__pycache__" EXCLUDE)
    install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/rsml DESTINATION ${CPLANTBOX_PYTHON_INSTALL_DIR}
            FILES_MATCHING
            PATTERN "*.py"
            PATTERN "__pycache__" EXCLUDE)
    install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/structural DESTINATION ${CPLANTBOX_PYTHON_INSTALL_DIR}
            FILES_MATCHING
            PATTERN "*.py"
            PATTERN "__pycache__" EXCLUDE)
    install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/functional DESTINATION ${CPLANTBOX_PYTHON_INSTALL_DIR}
            FILES_MATCHING
            PATTERN "*.py"
            PATTERN "__pycache__" EXCLUDE)
    install(DIRECTORY ${PROJECT_SOURCE_DIR}/modelparameter DESTINATION ${CPLANTBOX_PYTHON_INSTALL_DIR}
            PATTERN "__pycache__" EXCLUDE
            PATTERN "*.ipynb" EXCLUDE
            PATTERN ".DS_Store" EXCLUDE)
endif()
