# Build the SuiteSparse static libraries CPlantBox needs (SuiteSparse_config,
# AMD, BTF, COLAMD, KLU) with MSVC, mirroring the upstream Lib/Makefiles.
#
# SuiteSparse 5.x only ships Makefiles, which MSVC cannot use. The Makefiles
# compile most sources several times with different preprocessor definitions
# (integer and long-integer variants, real and complex variants). This project
# reproduces that exactly by generating tiny stub sources that set the
# definitions and then include the pinned upstream source file.
#
# Usage:
#   cmake -S scripts/deps/suitesparse-msvc -B <build> \
#     -DSUITESPARSE_SOURCE_DIR=<checkout of SuiteSparse v5.3.0>
#   cmake --build <build> --config Release

cmake_minimum_required(VERSION 3.18)
project(cplantbox-suitesparse-msvc C)

if(NOT DEFINED SUITESPARSE_SOURCE_DIR)
  message(FATAL_ERROR "Pass -DSUITESPARSE_SOURCE_DIR=<SuiteSparse checkout>")
endif()
get_filename_component(SS "${SUITESPARSE_SOURCE_DIR}" ABSOLUTE)
if(NOT EXISTS "${SS}/SuiteSparse_config/SuiteSparse_config.h")
  message(FATAL_ERROR "${SS} does not look like a SuiteSparse checkout")
endif()

set(STUB_DIR "${CMAKE_BINARY_DIR}/stubs")
file(MAKE_DIRECTORY "${STUB_DIR}")

# NTIMER: SuiteSparse_config's POSIX timer is not available under MSVC and is
# unused by CPlantBox; the upstream sources support opting out with -DNTIMER.
add_compile_definitions(NTIMER)

# Write a stub that applies `defs` and includes the real source.
# _ss_stub(<out-var> <stub-name> <relative-source-path> [defs...])
function(_ss_stub out_var stub_name source_rel)
  set(_stub "${STUB_DIR}/${stub_name}.c")
  set(_content "")
  foreach(_def ${ARGN})
    string(APPEND _content "#define ${_def}\n")
  endforeach()
  string(APPEND _content "#include \"${SS}/${source_rel}\"\n")
  file(WRITE "${_stub}" "${_content}")
  set(${out_var} "${_stub}" PARENT_SCOPE)
endfunction()

# ---------------------------------------------------------------------------
# SuiteSparse_config: single variant
# ---------------------------------------------------------------------------
add_library(suitesparseconfig STATIC "${SS}/SuiteSparse_config/SuiteSparse_config.c")
target_include_directories(suitesparseconfig PUBLIC "${SS}/SuiteSparse_config")

# ---------------------------------------------------------------------------
# AMD: every source compiled with -DDINT and again with -DDLONG
# (AMD/Lib/Makefile: amd_i_%.o / amd_l_%.o)
# ---------------------------------------------------------------------------
set(AMD_SOURCES amd_aat amd_1 amd_2 amd_dump amd_postorder amd_defaults
    amd_post_tree amd_order amd_control amd_info amd_valid amd_preprocess)
set(_amd_stubs "")
foreach(_f ${AMD_SOURCES})
  _ss_stub(_s "amd_i_${_f}" "AMD/Source/${_f}.c" DINT)
  list(APPEND _amd_stubs "${_s}")
  _ss_stub(_s "amd_l_${_f}" "AMD/Source/${_f}.c" DLONG)
  list(APPEND _amd_stubs "${_s}")
endforeach()
add_library(amd STATIC ${_amd_stubs})
target_include_directories(amd PRIVATE "${SS}/AMD/Include" "${SS}/SuiteSparse_config")

# ---------------------------------------------------------------------------
# BTF: every source compiled plain (int) and with -DDLONG
# (BTF/Lib/Makefile: btf_*.o / btf_l_*.o)
# ---------------------------------------------------------------------------
set(BTF_SOURCES btf_order btf_maxtrans btf_strongcomp)
set(_btf_stubs "")
foreach(_f ${BTF_SOURCES})
  _ss_stub(_s "btf_i_${_f}" "BTF/Source/${_f}.c")
  list(APPEND _btf_stubs "${_s}")
  _ss_stub(_s "btf_l_${_f}" "BTF/Source/${_f}.c" DLONG)
  list(APPEND _btf_stubs "${_s}")
endforeach()
add_library(btf STATIC ${_btf_stubs})
target_include_directories(btf PRIVATE "${SS}/BTF/Include" "${SS}/SuiteSparse_config")

# ---------------------------------------------------------------------------
# COLAMD: colamd.c compiled plain and with -DDLONG
# (COLAMD/Lib/Makefile: colamd.o / colamd_l.o)
# ---------------------------------------------------------------------------
_ss_stub(_colamd_i "colamd_i" "COLAMD/Source/colamd.c")
_ss_stub(_colamd_l "colamd_l" "COLAMD/Source/colamd.c" DLONG)
add_library(colamd STATIC "${_colamd_i}" "${_colamd_l}")
target_include_directories(colamd PRIVATE "${SS}/COLAMD/Include" "${SS}/SuiteSparse_config")

# ---------------------------------------------------------------------------
# KLU
# (KLU/Lib/Makefile: core sources in four variants {d, z, l, zl}; the common
# sources in two variants {int, long})
# ---------------------------------------------------------------------------
set(KLU_CORE_SOURCES klu klu_kernel klu_dump klu_factor klu_free_numeric
    klu_solve klu_scale klu_refactor klu_tsolve klu_diagnostics klu_sort
    klu_extract)
set(KLU_COMMON_SOURCES klu_free_symbolic klu_defaults klu_analyze_given
    klu_analyze klu_memory)
set(_klu_stubs "")
foreach(_f ${KLU_CORE_SOURCES})
  _ss_stub(_s "${_f}_d" "KLU/Source/${_f}.c")
  list(APPEND _klu_stubs "${_s}")
  _ss_stub(_s "${_f}_z" "KLU/Source/${_f}.c" COMPLEX)
  list(APPEND _klu_stubs "${_s}")
  _ss_stub(_s "${_f}_l" "KLU/Source/${_f}.c" DLONG)
  list(APPEND _klu_stubs "${_s}")
  _ss_stub(_s "${_f}_zl" "KLU/Source/${_f}.c" COMPLEX DLONG)
  list(APPEND _klu_stubs "${_s}")
endforeach()
foreach(_f ${KLU_COMMON_SOURCES})
  _ss_stub(_s "${_f}_i" "KLU/Source/${_f}.c")
  list(APPEND _klu_stubs "${_s}")
  _ss_stub(_s "${_f}_l" "KLU/Source/${_f}.c" DLONG)
  list(APPEND _klu_stubs "${_s}")
endforeach()
add_library(klu STATIC ${_klu_stubs})
target_include_directories(klu PRIVATE
  "${SS}/KLU/Include" "${SS}/AMD/Include" "${SS}/COLAMD/Include"
  "${SS}/BTF/Include" "${SS}/SuiteSparse_config")
