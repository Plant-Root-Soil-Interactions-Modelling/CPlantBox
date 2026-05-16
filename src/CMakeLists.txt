# -----------------------------
# Build options
# -----------------------------
option(CPB_BUILD_PYTHON "Build the Python bindings" ON)
option(CPB_BUILD_SHARED "Build CPlantBox as a shared library" ON)

message(STATUS "Project dir: ${PROJECT_SOURCE_DIR}")
message(STATUS "Current source dir: ${CMAKE_CURRENT_SOURCE_DIR}")
message(STATUS "Project binary dir: ${PROJECT_BINARY_DIR}")

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
    INSTALL_RPATH "$ORIGIN"
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
function(cpb_add_imported_archive target archive_path)
    add_library(${target} STATIC IMPORTED)

    set(_binary_archive "${CMAKE_CURRENT_BINARY_DIR}/${archive_path}")
    set(_source_archive "${CMAKE_CURRENT_SOURCE_DIR}/${archive_path}")

    if(EXISTS "${_binary_archive}")
        set(_archive_location "${_binary_archive}")
    elseif(EXISTS "${_source_archive}")
        set(_archive_location "${_source_archive}")
    else()
        message(FATAL_ERROR "Missing native dependency archive for ${target}: ${_binary_archive} or ${_source_archive}")
    endif()

    set_target_properties(${target} PROPERTIES IMPORTED_LOCATION "${_archive_location}")
endfunction()

if(PIAFMUNCH)
    target_compile_definitions(CPlantBox PUBLIC ENABLE_PIAFMUNCH)

    cpb_add_imported_archive(libklu external/suitsparse/lib/libklu.a)
    cpb_add_imported_archive(libamd external/suitsparse/lib/libamd.a)
    cpb_add_imported_archive(libbtf external/suitsparse/lib/libbtf.a)
    cpb_add_imported_archive(libcolamd external/suitsparse/lib/libcolamd.a)
    cpb_add_imported_archive(suitesparseconfig external/suitsparse/lib/libsuitesparseconfig.a)

    cpb_add_imported_archive(arkode external/sundials/lib/libsundials_arkode.a)
    cpb_add_imported_archive(sundials_cvode external/sundials/lib/libsundials_cvode.a)
    cpb_add_imported_archive(sundials_nvecserial external/sundials/lib/libsundials_nvecserial.a)
    cpb_add_imported_archive(sundials_sunlinsolband external/sundials/lib/libsundials_sunlinsolband.a)
    cpb_add_imported_archive(sundials_sunlinsoldense external/sundials/lib/libsundials_sunlinsoldense.a)
    cpb_add_imported_archive(sundials_sunlinsolklu external/sundials/lib/libsundials_sunlinsolklu.a)
    cpb_add_imported_archive(sundials_sunlinsolpcg external/sundials/lib/libsundials_sunlinsolpcg.a)
    cpb_add_imported_archive(sundials_sunlinsolspbcgs external/sundials/lib/libsundials_sunlinsolspbcgs.a)
    cpb_add_imported_archive(sundials_sunlinsolspfgmr external/sundials/lib/libsundials_sunlinsolspfgmr.a)
    cpb_add_imported_archive(sundials_sunlinsolspgmr external/sundials/lib/libsundials_sunlinsolspgmr.a)
    cpb_add_imported_archive(sundials_sunlinsolsptfqmr external/sundials/lib/libsundials_sunlinsolsptfqmr.a)
    cpb_add_imported_archive(sundials_sunmatrixband external/sundials/lib/libsundials_sunmatrixband.a)
    cpb_add_imported_archive(sundials_sunmatrixdense external/sundials/lib/libsundials_sunmatrixdense.a)
    cpb_add_imported_archive(sundials_sunmatrixsparse external/sundials/lib/libsundials_sunmatrixsparse.a)
    cpb_add_imported_archive(sundials_sunnonlinsolnewton external/sundials/lib/libsundials_sunnonlinsolnewton.a)

    # Static archive order matters for one-pass linkers: users first, then the
    # archives that provide their unresolved symbols.
    target_link_libraries(CPlantBox
        PRIVATE
            arkode
            sundials_cvode
            sundials_sunlinsolklu
            sundials_sunlinsolspgmr
            sundials_sunlinsolspfgmr
            sundials_sunlinsolspbcgs
            sundials_sunlinsolpcg
            sundials_sunlinsolsptfqmr
            sundials_sunlinsolband
            sundials_sunlinsoldense
            sundials_sunnonlinsolnewton
            sundials_sunmatrixsparse
            sundials_sunmatrixdense
            sundials_sunmatrixband
            sundials_nvecserial
            libklu
            libbtf
            libcolamd
            libamd
            suitesparseconfig
    )

    target_include_directories(CPlantBox
        PUBLIC
            ${CMAKE_CURRENT_SOURCE_DIR}/external/sundials/include
            ${CMAKE_CURRENT_SOURCE_DIR}/external/suitsparse/include
    )
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
        INSTALL_RPATH "$ORIGIN"
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
