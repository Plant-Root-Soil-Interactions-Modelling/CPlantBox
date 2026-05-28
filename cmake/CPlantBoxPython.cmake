# CPlantBox Python extension and package install configuration.
#
# Keep pybind11 module creation and Python package install rules outside the
# core native target definition so src/CMakeLists.txt remains target-oriented.

set(_CPB_PYTHON_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
set(_CPB_PYTHON_BINARY_DIR "${CMAKE_CURRENT_BINARY_DIR}")
set(_CPB_PYTHON_PROJECT_SOURCE_DIR "${PROJECT_SOURCE_DIR}")

function(cplantbox_configure_python_bindings target)
    if(NOT TARGET "${target}")
        message(FATAL_ERROR "cplantbox_configure_python_bindings called with unknown target: ${target}")
    endif()

    if(CMAKE_VERSION VERSION_LESS 3.18)
        find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
    else()
        find_package(Python3 COMPONENTS Interpreter Development.Module REQUIRED)
    endif()
    message(STATUS "Python library path (site-packages): ${Python3_SITEARCH}")

    add_subdirectory("${_CPB_PYTHON_SOURCE_DIR}/external/pybind11" "${_CPB_PYTHON_BINARY_DIR}/external/pybind11")

    pybind11_add_module(plantbox "${_CPB_PYTHON_SOURCE_DIR}/PyPlantBox.cpp")
    target_link_libraries(plantbox PRIVATE ${target})
    if(TARGET Python3::Module)
        target_link_libraries(plantbox PRIVATE Python3::Module)
    endif()

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
        install(TARGETS ${target}
                LIBRARY DESTINATION "${CPLANTBOX_PYTHON_INSTALL_DIR}"
                RUNTIME DESTINATION "${CPLANTBOX_PYTHON_INSTALL_DIR}")
    endif()
    install(TARGETS plantbox
            LIBRARY DESTINATION "${CPLANTBOX_PYTHON_INSTALL_DIR}"
            RUNTIME DESTINATION "${CPLANTBOX_PYTHON_INSTALL_DIR}")

    install(FILES ${_CPB_PYTHON_SOURCE_DIR}/__init__.py DESTINATION ${CPLANTBOX_PYTHON_INSTALL_DIR})

    install(DIRECTORY ${_CPB_PYTHON_SOURCE_DIR}/visualisation DESTINATION ${CPLANTBOX_PYTHON_INSTALL_DIR}
            FILES_MATCHING
            PATTERN "*.py"
            PATTERN "__pycache__" EXCLUDE)
    install(DIRECTORY ${_CPB_PYTHON_SOURCE_DIR}/rsml DESTINATION ${CPLANTBOX_PYTHON_INSTALL_DIR}
            FILES_MATCHING
            PATTERN "*.py"
            PATTERN "__pycache__" EXCLUDE)
    install(DIRECTORY ${_CPB_PYTHON_SOURCE_DIR}/structural DESTINATION ${CPLANTBOX_PYTHON_INSTALL_DIR}
            FILES_MATCHING
            PATTERN "*.py"
            PATTERN "__pycache__" EXCLUDE)
    install(DIRECTORY ${_CPB_PYTHON_SOURCE_DIR}/functional DESTINATION ${CPLANTBOX_PYTHON_INSTALL_DIR}
            FILES_MATCHING
            PATTERN "*.py"
            PATTERN "__pycache__" EXCLUDE)
    install(DIRECTORY ${_CPB_PYTHON_PROJECT_SOURCE_DIR}/modelparameter DESTINATION ${CPLANTBOX_PYTHON_INSTALL_DIR}
            PATTERN "__pycache__" EXCLUDE
            PATTERN "*.ipynb" EXCLUDE
            PATTERN ".DS_Store" EXCLUDE)

    if(SKBUILD AND CPB_NATIVE_DEPS_PREFIX AND SKBUILD_METADATA_DIR)
        set(CPLANTBOX_THIRD_PARTY_LICENSE_DIR
            "${SKBUILD_METADATA_DIR}/licenses/third_party/native-dependencies")
        install(DIRECTORY "${CPB_NATIVE_DEPS_PREFIX}/share/licenses/"
                DESTINATION "${CPLANTBOX_THIRD_PARTY_LICENSE_DIR}")
        install(FILES "${CPB_NATIVE_DEPS_PREFIX}/share/cplantbox-native-deps-provenance.txt"
                DESTINATION "${CPLANTBOX_THIRD_PARTY_LICENSE_DIR}")
    endif()
endfunction()
