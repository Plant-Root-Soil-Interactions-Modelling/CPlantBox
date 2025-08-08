cmake_minimum_required (VERSION 2.8.12)
project (CPlantBox)

# set compiler options
set(RELEASE_WARNING_OPTS "-Wall -Wunused -Wmissing-include-dirs -Wcast-align -Wno-sign-compare  ")
set(DEBUG_WARNING_OPTS "-Wall -Wunused -Wmissing-include-dirs -Wcast-align -Wno-sign-compare ")
# Portable release opts (no -march=native for wheel portability)
set(RELEASE_OPTS "-std=c++17 -fno-strict-aliasing -fstrict-overflow -fno-finite-math-only -O3 -funroll-loops -g0 ")
set(CMAKE_CXX_FLAGS_RELEASE "${RELEASE_WARNING_OPTS} ${RELEASE_OPTS}")
set(CMAKE_CXX_FLAGS_DEBUG "-std=c++17 -O0 -ggdb -Wall")

# Default build type: Release (can be overridden by user/toolchain)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# add source directory to the include path
include_directories(${PROJECT_SOURCE_DIR}/src/structural)
include_directories(${PROJECT_SOURCE_DIR}/src/functional)

# add subdirectories
add_subdirectory(src)
#add_subdirectory(tutorial)

# build tutorial (TODO)
# INCLUDE(UseLATEX.cmake)
