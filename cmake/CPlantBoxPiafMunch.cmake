# CPlantBox PiafMunch/native dependency configuration.
#
# This module intentionally keeps the SuiteSparse/SUNDIALS provider policy and
# static link order in one place so src/CMakeLists.txt can stay target-oriented.

set(_CPB_PIAFMUNCH_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
set(_CPB_PIAFMUNCH_BINARY_DIR "${CMAKE_CURRENT_BINARY_DIR}")

function(cpb_validate_provider variable)
    if(NOT "${${variable}}" STREQUAL "bundled" AND NOT "${${variable}}" STREQUAL "source" AND NOT "${${variable}}" STREQUAL "system")
        message(FATAL_ERROR "${variable} must be 'bundled', 'source', or 'system', got '${${variable}}'")
    endif()
endfunction()

function(cpb_add_imported_library target library_type library_location)
    add_library(${target} ${library_type} IMPORTED)
    set_target_properties(${target} PROPERTIES IMPORTED_LOCATION "${library_location}")
endfunction()

function(cpb_add_bundled_archive target archive_path)
    set(_binary_archive "${_CPB_PIAFMUNCH_BINARY_DIR}/${archive_path}")
    set(_source_archive "${_CPB_PIAFMUNCH_SOURCE_DIR}/${archive_path}")

    if(EXISTS "${_binary_archive}")
        set(_archive_location "${_binary_archive}")
    elseif(EXISTS "${_source_archive}")
        set(_archive_location "${_source_archive}")
    else()
        message(FATAL_ERROR "Missing bundled native dependency archive for ${target}: ${_binary_archive} or ${_source_archive}")
    endif()

    cpb_add_imported_library(${target} STATIC "${_archive_location}")
endfunction()

function(cpb_find_system_library output_variable display_name)
    find_library(_library_path NAMES ${ARGN})
    if(NOT _library_path)
        message(FATAL_ERROR "Requested system ${display_name}, but none of these libraries were found: ${ARGN}")
    endif()
    set(${output_variable} "${_library_path}" PARENT_SCOPE)
    unset(_library_path CACHE)
endfunction()

function(cpb_source_prefixes output_variable)
    set(_prefixes ${CPB_NATIVE_DEPS_PREFIX} ${CMAKE_PREFIX_PATH})
    list(REMOVE_ITEM _prefixes "")
    if(NOT _prefixes)
        message(FATAL_ERROR "Requested source-built native dependencies, but CPB_NATIVE_DEPS_PREFIX/CMAKE_PREFIX_PATH is empty")
    endif()
    set(${output_variable} ${_prefixes} PARENT_SCOPE)
endfunction()

function(cpb_find_source_header_dir output_variable display_name header)
    cpb_source_prefixes(_prefixes)
    find_path(_include_dir NAMES ${header} PATHS ${_prefixes} PATH_SUFFIXES include NO_DEFAULT_PATH)
    if(NOT _include_dir)
        message(FATAL_ERROR "Requested source-built ${display_name}, but ${header} was not found under: ${_prefixes}")
    endif()
    set(${output_variable} "${_include_dir}" PARENT_SCOPE)
    unset(_include_dir CACHE)
endfunction()

function(cpb_find_source_library output_variable display_name)
    cpb_source_prefixes(_prefixes)
    find_library(_library_path NAMES ${ARGN} PATHS ${_prefixes} PATH_SUFFIXES lib NO_DEFAULT_PATH)
    if(NOT _library_path)
        message(FATAL_ERROR "Requested source-built ${display_name}, but none of these libraries were found under ${_prefixes}: ${ARGN}")
    endif()
    if(NOT _library_path MATCHES "\\.a$")
        message(FATAL_ERROR "Requested source-built ${display_name}, but expected a static .a archive and found: ${_library_path}")
    endif()
    set(${output_variable} "${_library_path}" PARENT_SCOPE)
    unset(_library_path CACHE)
endfunction()

function(cpb_require_header include_dir header display_name)
    if(NOT EXISTS "${include_dir}/${header}")
        message(FATAL_ERROR "Missing ${display_name} header: ${include_dir}/${header}")
    endif()
endfunction()

function(cplantbox_configure_piafmunch target)
    if(NOT TARGET "${target}")
        message(FATAL_ERROR "cplantbox_configure_piafmunch called with unknown target: ${target}")
    endif()

    cpb_validate_provider(CPB_SUITESPARSE_PROVIDER)
    cpb_validate_provider(CPB_SUNDIALS_PROVIDER)

    message(STATUS "CPlantBox SuiteSparse provider: ${CPB_SUITESPARSE_PROVIDER}")
    message(STATUS "CPlantBox SUNDIALS provider: ${CPB_SUNDIALS_PROVIDER}")

    target_compile_definitions(${target} PUBLIC ENABLE_PIAFMUNCH)

    if(CPB_SUITESPARSE_PROVIDER STREQUAL "bundled")
        set(CPB_SUITESPARSE_INCLUDE_DIR "${_CPB_PIAFMUNCH_SOURCE_DIR}/external/suitsparse/include")
        cpb_add_bundled_archive(libklu external/suitsparse/lib/libklu.a)
        cpb_add_bundled_archive(libamd external/suitsparse/lib/libamd.a)
        cpb_add_bundled_archive(libbtf external/suitsparse/lib/libbtf.a)
        cpb_add_bundled_archive(libcolamd external/suitsparse/lib/libcolamd.a)
        cpb_add_bundled_archive(suitesparseconfig external/suitsparse/lib/libsuitesparseconfig.a)
    elseif(CPB_SUITESPARSE_PROVIDER STREQUAL "source")
        cpb_find_source_header_dir(CPB_SUITESPARSE_INCLUDE_DIR "SuiteSparse" klu.h)
        cpb_find_source_library(CPB_LIBKLU "SuiteSparse KLU" klu libklu)
        cpb_find_source_library(CPB_LIBAMD "SuiteSparse AMD" amd libamd)
        cpb_find_source_library(CPB_LIBBTF "SuiteSparse BTF" btf libbtf)
        cpb_find_source_library(CPB_LIBCOLAMD "SuiteSparse COLAMD" colamd libcolamd)
        cpb_find_source_library(CPB_LIBSUITESPARSECONFIG "SuiteSparse config" suitesparseconfig libsuitesparseconfig)

        cpb_require_header("${CPB_SUITESPARSE_INCLUDE_DIR}" klu.h "SuiteSparse")
        cpb_require_header("${CPB_SUITESPARSE_INCLUDE_DIR}" amd.h "SuiteSparse")
        cpb_require_header("${CPB_SUITESPARSE_INCLUDE_DIR}" btf.h "SuiteSparse")
        cpb_require_header("${CPB_SUITESPARSE_INCLUDE_DIR}" colamd.h "SuiteSparse")
        cpb_require_header("${CPB_SUITESPARSE_INCLUDE_DIR}" SuiteSparse_config.h "SuiteSparse")

        cpb_add_imported_library(libklu STATIC "${CPB_LIBKLU}")
        cpb_add_imported_library(libamd STATIC "${CPB_LIBAMD}")
        cpb_add_imported_library(libbtf STATIC "${CPB_LIBBTF}")
        cpb_add_imported_library(libcolamd STATIC "${CPB_LIBCOLAMD}")
        cpb_add_imported_library(suitesparseconfig STATIC "${CPB_LIBSUITESPARSECONFIG}")
    else()
        find_path(CPB_SUITESPARSE_INCLUDE_DIR NAMES klu.h PATH_SUFFIXES suitesparse)
        if(NOT CPB_SUITESPARSE_INCLUDE_DIR)
            message(FATAL_ERROR "Requested system SuiteSparse, but required headers were not found")
        endif()

        cpb_find_system_library(CPB_LIBKLU "SuiteSparse KLU" klu libklu)
        cpb_find_system_library(CPB_LIBAMD "SuiteSparse AMD" amd libamd)
        cpb_find_system_library(CPB_LIBBTF "SuiteSparse BTF" btf libbtf)
        cpb_find_system_library(CPB_LIBCOLAMD "SuiteSparse COLAMD" colamd libcolamd)
        cpb_find_system_library(CPB_LIBSUITESPARSECONFIG "SuiteSparse config" suitesparseconfig libsuitesparseconfig)

        cpb_require_header("${CPB_SUITESPARSE_INCLUDE_DIR}" klu.h "SuiteSparse")
        cpb_require_header("${CPB_SUITESPARSE_INCLUDE_DIR}" amd.h "SuiteSparse")
        cpb_require_header("${CPB_SUITESPARSE_INCLUDE_DIR}" btf.h "SuiteSparse")
        cpb_require_header("${CPB_SUITESPARSE_INCLUDE_DIR}" colamd.h "SuiteSparse")
        cpb_require_header("${CPB_SUITESPARSE_INCLUDE_DIR}" SuiteSparse_config.h "SuiteSparse")

        cpb_add_imported_library(libklu UNKNOWN "${CPB_LIBKLU}")
        cpb_add_imported_library(libamd UNKNOWN "${CPB_LIBAMD}")
        cpb_add_imported_library(libbtf UNKNOWN "${CPB_LIBBTF}")
        cpb_add_imported_library(libcolamd UNKNOWN "${CPB_LIBCOLAMD}")
        cpb_add_imported_library(suitesparseconfig UNKNOWN "${CPB_LIBSUITESPARSECONFIG}")
    endif()

    if(CPB_SUNDIALS_PROVIDER STREQUAL "bundled")
        set(CPB_SUNDIALS_INCLUDE_DIR "${_CPB_PIAFMUNCH_SOURCE_DIR}/external/sundials/include")
        cpb_add_bundled_archive(arkode external/sundials/lib/libsundials_arkode.a)
        cpb_add_bundled_archive(sundials_cvode external/sundials/lib/libsundials_cvode.a)
        cpb_add_bundled_archive(sundials_nvecserial external/sundials/lib/libsundials_nvecserial.a)
        cpb_add_bundled_archive(sundials_sunlinsolband external/sundials/lib/libsundials_sunlinsolband.a)
        cpb_add_bundled_archive(sundials_sunlinsoldense external/sundials/lib/libsundials_sunlinsoldense.a)
        cpb_add_bundled_archive(sundials_sunlinsolklu external/sundials/lib/libsundials_sunlinsolklu.a)
        cpb_add_bundled_archive(sundials_sunlinsolpcg external/sundials/lib/libsundials_sunlinsolpcg.a)
        cpb_add_bundled_archive(sundials_sunlinsolspbcgs external/sundials/lib/libsundials_sunlinsolspbcgs.a)
        cpb_add_bundled_archive(sundials_sunlinsolspfgmr external/sundials/lib/libsundials_sunlinsolspfgmr.a)
        cpb_add_bundled_archive(sundials_sunlinsolspgmr external/sundials/lib/libsundials_sunlinsolspgmr.a)
        cpb_add_bundled_archive(sundials_sunlinsolsptfqmr external/sundials/lib/libsundials_sunlinsolsptfqmr.a)
        cpb_add_bundled_archive(sundials_sunmatrixband external/sundials/lib/libsundials_sunmatrixband.a)
        cpb_add_bundled_archive(sundials_sunmatrixdense external/sundials/lib/libsundials_sunmatrixdense.a)
        cpb_add_bundled_archive(sundials_sunmatrixsparse external/sundials/lib/libsundials_sunmatrixsparse.a)
        cpb_add_bundled_archive(sundials_sunnonlinsolnewton external/sundials/lib/libsundials_sunnonlinsolnewton.a)
    elseif(CPB_SUNDIALS_PROVIDER STREQUAL "source")
        cpb_find_source_header_dir(CPB_SUNDIALS_INCLUDE_DIR "SUNDIALS" sundials/sundials_config.h)
        cpb_find_source_library(CPB_LIBSUNDIALS_ARKODE "SUNDIALS ARKODE" sundials_arkode libsundials_arkode)
        cpb_find_source_library(CPB_LIBSUNDIALS_CVODE "SUNDIALS CVODE" sundials_cvode libsundials_cvode)
        cpb_find_source_library(CPB_LIBSUNDIALS_NVEC "SUNDIALS serial nvector" sundials_nvecserial libsundials_nvecserial)
        cpb_find_source_library(CPB_LIBSUNDIALS_SUNLINSOLBAND "SUNDIALS band linear solver" sundials_sunlinsolband libsundials_sunlinsolband)
        cpb_find_source_library(CPB_LIBSUNDIALS_SUNLINSOLDENSE "SUNDIALS dense linear solver" sundials_sunlinsoldense libsundials_sunlinsoldense)
        cpb_find_source_library(CPB_LIBSUNDIALS_SUNLINSOLKLU "SUNDIALS KLU linear solver" sundials_sunlinsolklu libsundials_sunlinsolklu)
        cpb_find_source_library(CPB_LIBSUNDIALS_SUNLINSOLPCG "SUNDIALS PCG linear solver" sundials_sunlinsolpcg libsundials_sunlinsolpcg)
        cpb_find_source_library(CPB_LIBSUNDIALS_SUNLINSOLSPBCGS "SUNDIALS SPBCGS linear solver" sundials_sunlinsolspbcgs libsundials_sunlinsolspbcgs)
        cpb_find_source_library(CPB_LIBSUNDIALS_SUNLINSOLSPFGMR "SUNDIALS SPFGMR linear solver" sundials_sunlinsolspfgmr libsundials_sunlinsolspfgmr)
        cpb_find_source_library(CPB_LIBSUNDIALS_SUNLINSOLSPGMR "SUNDIALS SPGMR linear solver" sundials_sunlinsolspgmr libsundials_sunlinsolspgmr)
        cpb_find_source_library(CPB_LIBSUNDIALS_SUNLINSOLSPTFQMR "SUNDIALS SPTFQMR linear solver" sundials_sunlinsolsptfqmr libsundials_sunlinsolsptfqmr)
        cpb_find_source_library(CPB_LIBSUNDIALS_SUNMATRIXBAND "SUNDIALS band matrix" sundials_sunmatrixband libsundials_sunmatrixband)
        cpb_find_source_library(CPB_LIBSUNDIALS_SUNMATRIXDENSE "SUNDIALS dense matrix" sundials_sunmatrixdense libsundials_sunmatrixdense)
        cpb_find_source_library(CPB_LIBSUNDIALS_SUNMATRIXSPARSE "SUNDIALS sparse matrix" sundials_sunmatrixsparse libsundials_sunmatrixsparse)
        cpb_find_source_library(CPB_LIBSUNDIALS_SUNNONLINSOLNEWTON "SUNDIALS Newton nonlinear solver" sundials_sunnonlinsolnewton libsundials_sunnonlinsolnewton)

        cpb_require_header("${CPB_SUNDIALS_INCLUDE_DIR}" sundials/sundials_config.h "SUNDIALS")
        cpb_require_header("${CPB_SUNDIALS_INCLUDE_DIR}" cvode/cvode.h "SUNDIALS")
        cpb_require_header("${CPB_SUNDIALS_INCLUDE_DIR}" nvector/nvector_serial.h "SUNDIALS")
        cpb_require_header("${CPB_SUNDIALS_INCLUDE_DIR}" sunlinsol/sunlinsol_klu.h "SUNDIALS")
        cpb_require_header("${CPB_SUNDIALS_INCLUDE_DIR}" sunmatrix/sunmatrix_sparse.h "SUNDIALS")

        cpb_add_imported_library(arkode STATIC "${CPB_LIBSUNDIALS_ARKODE}")
        cpb_add_imported_library(sundials_cvode STATIC "${CPB_LIBSUNDIALS_CVODE}")
        cpb_add_imported_library(sundials_nvecserial STATIC "${CPB_LIBSUNDIALS_NVEC}")
        cpb_add_imported_library(sundials_sunlinsolband STATIC "${CPB_LIBSUNDIALS_SUNLINSOLBAND}")
        cpb_add_imported_library(sundials_sunlinsoldense STATIC "${CPB_LIBSUNDIALS_SUNLINSOLDENSE}")
        cpb_add_imported_library(sundials_sunlinsolklu STATIC "${CPB_LIBSUNDIALS_SUNLINSOLKLU}")
        cpb_add_imported_library(sundials_sunlinsolpcg STATIC "${CPB_LIBSUNDIALS_SUNLINSOLPCG}")
        cpb_add_imported_library(sundials_sunlinsolspbcgs STATIC "${CPB_LIBSUNDIALS_SUNLINSOLSPBCGS}")
        cpb_add_imported_library(sundials_sunlinsolspfgmr STATIC "${CPB_LIBSUNDIALS_SUNLINSOLSPFGMR}")
        cpb_add_imported_library(sundials_sunlinsolspgmr STATIC "${CPB_LIBSUNDIALS_SUNLINSOLSPGMR}")
        cpb_add_imported_library(sundials_sunlinsolsptfqmr STATIC "${CPB_LIBSUNDIALS_SUNLINSOLSPTFQMR}")
        cpb_add_imported_library(sundials_sunmatrixband STATIC "${CPB_LIBSUNDIALS_SUNMATRIXBAND}")
        cpb_add_imported_library(sundials_sunmatrixdense STATIC "${CPB_LIBSUNDIALS_SUNMATRIXDENSE}")
        cpb_add_imported_library(sundials_sunmatrixsparse STATIC "${CPB_LIBSUNDIALS_SUNMATRIXSPARSE}")
        cpb_add_imported_library(sundials_sunnonlinsolnewton STATIC "${CPB_LIBSUNDIALS_SUNNONLINSOLNEWTON}")
    else()
        find_path(CPB_SUNDIALS_INCLUDE_DIR NAMES sundials/sundials_config.h)
        if(NOT CPB_SUNDIALS_INCLUDE_DIR)
            message(FATAL_ERROR "Requested system SUNDIALS, but required headers were not found")
        endif()

        cpb_find_system_library(CPB_LIBSUNDIALS_ARKODE "SUNDIALS ARKODE" sundials_arkode libsundials_arkode)
        cpb_find_system_library(CPB_LIBSUNDIALS_CVODE "SUNDIALS CVODE" sundials_cvode libsundials_cvode)
        cpb_find_system_library(CPB_LIBSUNDIALS_NVEC "SUNDIALS serial nvector" sundials_nvecserial libsundials_nvecserial)
        cpb_find_system_library(CPB_LIBSUNDIALS_SUNLINSOLBAND "SUNDIALS band linear solver" sundials_sunlinsolband libsundials_sunlinsolband)
        cpb_find_system_library(CPB_LIBSUNDIALS_SUNLINSOLDENSE "SUNDIALS dense linear solver" sundials_sunlinsoldense libsundials_sunlinsoldense)
        cpb_find_system_library(CPB_LIBSUNDIALS_SUNLINSOLKLU "SUNDIALS KLU linear solver" sundials_sunlinsolklu libsundials_sunlinsolklu)
        cpb_find_system_library(CPB_LIBSUNDIALS_SUNLINSOLPCG "SUNDIALS PCG linear solver" sundials_sunlinsolpcg libsundials_sunlinsolpcg)
        cpb_find_system_library(CPB_LIBSUNDIALS_SUNLINSOLSPBCGS "SUNDIALS SPBCGS linear solver" sundials_sunlinsolspbcgs libsundials_sunlinsolspbcgs)
        cpb_find_system_library(CPB_LIBSUNDIALS_SUNLINSOLSPFGMR "SUNDIALS SPFGMR linear solver" sundials_sunlinsolspfgmr libsundials_sunlinsolspfgmr)
        cpb_find_system_library(CPB_LIBSUNDIALS_SUNLINSOLSPGMR "SUNDIALS SPGMR linear solver" sundials_sunlinsolspgmr libsundials_sunlinsolspgmr)
        cpb_find_system_library(CPB_LIBSUNDIALS_SUNLINSOLSPTFQMR "SUNDIALS SPTFQMR linear solver" sundials_sunlinsolsptfqmr libsundials_sunlinsolsptfqmr)
        cpb_find_system_library(CPB_LIBSUNDIALS_SUNMATRIXBAND "SUNDIALS band matrix" sundials_sunmatrixband libsundials_sunmatrixband)
        cpb_find_system_library(CPB_LIBSUNDIALS_SUNMATRIXDENSE "SUNDIALS dense matrix" sundials_sunmatrixdense libsundials_sunmatrixdense)
        cpb_find_system_library(CPB_LIBSUNDIALS_SUNMATRIXSPARSE "SUNDIALS sparse matrix" sundials_sunmatrixsparse libsundials_sunmatrixsparse)
        cpb_find_system_library(CPB_LIBSUNDIALS_SUNNONLINSOLNEWTON "SUNDIALS Newton nonlinear solver" sundials_sunnonlinsolnewton libsundials_sunnonlinsolnewton)

        cpb_require_header("${CPB_SUNDIALS_INCLUDE_DIR}" sundials/sundials_config.h "SUNDIALS")
        cpb_require_header("${CPB_SUNDIALS_INCLUDE_DIR}" cvode/cvode.h "SUNDIALS")
        cpb_require_header("${CPB_SUNDIALS_INCLUDE_DIR}" nvector/nvector_serial.h "SUNDIALS")
        cpb_require_header("${CPB_SUNDIALS_INCLUDE_DIR}" sunlinsol/sunlinsol_klu.h "SUNDIALS")
        cpb_require_header("${CPB_SUNDIALS_INCLUDE_DIR}" sunmatrix/sunmatrix_sparse.h "SUNDIALS")

        cpb_add_imported_library(arkode UNKNOWN "${CPB_LIBSUNDIALS_ARKODE}")
        cpb_add_imported_library(sundials_cvode UNKNOWN "${CPB_LIBSUNDIALS_CVODE}")
        cpb_add_imported_library(sundials_nvecserial UNKNOWN "${CPB_LIBSUNDIALS_NVEC}")
        cpb_add_imported_library(sundials_sunlinsolband UNKNOWN "${CPB_LIBSUNDIALS_SUNLINSOLBAND}")
        cpb_add_imported_library(sundials_sunlinsoldense UNKNOWN "${CPB_LIBSUNDIALS_SUNLINSOLDENSE}")
        cpb_add_imported_library(sundials_sunlinsolklu UNKNOWN "${CPB_LIBSUNDIALS_SUNLINSOLKLU}")
        cpb_add_imported_library(sundials_sunlinsolpcg UNKNOWN "${CPB_LIBSUNDIALS_SUNLINSOLPCG}")
        cpb_add_imported_library(sundials_sunlinsolspbcgs UNKNOWN "${CPB_LIBSUNDIALS_SUNLINSOLSPBCGS}")
        cpb_add_imported_library(sundials_sunlinsolspfgmr UNKNOWN "${CPB_LIBSUNDIALS_SUNLINSOLSPFGMR}")
        cpb_add_imported_library(sundials_sunlinsolspgmr UNKNOWN "${CPB_LIBSUNDIALS_SUNLINSOLSPGMR}")
        cpb_add_imported_library(sundials_sunlinsolsptfqmr UNKNOWN "${CPB_LIBSUNDIALS_SUNLINSOLSPTFQMR}")
        cpb_add_imported_library(sundials_sunmatrixband UNKNOWN "${CPB_LIBSUNDIALS_SUNMATRIXBAND}")
        cpb_add_imported_library(sundials_sunmatrixdense UNKNOWN "${CPB_LIBSUNDIALS_SUNMATRIXDENSE}")
        cpb_add_imported_library(sundials_sunmatrixsparse UNKNOWN "${CPB_LIBSUNDIALS_SUNMATRIXSPARSE}")
        cpb_add_imported_library(sundials_sunnonlinsolnewton UNKNOWN "${CPB_LIBSUNDIALS_SUNNONLINSOLNEWTON}")

        find_library(CPB_LIBSUNDIALS_CORE NAMES sundials_core libsundials_core)
        if(CPB_LIBSUNDIALS_CORE)
            cpb_add_imported_library(sundials_core UNKNOWN "${CPB_LIBSUNDIALS_CORE}")
        endif()
    endif()

    # Static archive order matters for one-pass linkers: users first, then the
    # archives that provide their unresolved symbols.
    target_link_libraries(${target}
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

    if(TARGET sundials_core)
        target_link_libraries(${target} PRIVATE sundials_core)
    endif()

    if(CPB_BUILD_SHARED)
        # Catch unresolved native symbols at link time instead of import time.
        target_link_libraries(${target} PRIVATE "$<$<PLATFORM_ID:Linux>:-Wl,-z,defs>")
    endif()

    target_include_directories(${target}
        PUBLIC
            ${CPB_SUNDIALS_INCLUDE_DIR}
            ${CPB_SUITESPARSE_INCLUDE_DIR}
    )
endfunction()
