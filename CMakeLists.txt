cmake_minimum_required (VERSION 3.12)
project (CPlantBox VERSION 2.1.0 LANGUAGES CXX)

option(PIAFMUNCH "Enable optional PIAFMUNCH code" ON)
if(SKBUILD)
    set(CPB_NATIVE_CPU_TUNING_DEFAULT OFF)
else()
    set(CPB_NATIVE_CPU_TUNING_DEFAULT ON)
endif()
option(CPB_NATIVE_CPU_TUNING "Enable native CPU tuning flags such as -march=native" ${CPB_NATIVE_CPU_TUNING_DEFAULT})

# set compiler options
set(RELEASE_WARNING_OPTS "-Wall -Wunused -Wmissing-include-dirs -Wcast-align -Wno-sign-compare  ")
set(DEBUG_WARNING_OPTS "-Wall -Wunused -Wmissing-include-dirs -Wcast-align -Wno-sign-compare ")
set(RELEASE_OPTS "-std=c++17 -fno-strict-aliasing -fstrict-overflow -fno-finite-math-only -O3 -funroll-loops -g0 ")
if(CPB_NATIVE_CPU_TUNING)
    set(RELEASE_OPTS "${RELEASE_OPTS} -march=native")
endif()
set(CMAKE_CXX_FLAGS_RELEASE "${RELEASE_WARNING_OPTS} ${RELEASE_OPTS}")
set(CMAKE_CXX_FLAGS_DEBUG "-std=c++17 -O0 -ggdb -Wall")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# set default build type to release Release, Debug, RelWithDebInfo
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Build type" FORCE)
endif()

# add subdirectories
add_subdirectory(src)
#add_subdirectory(tutorial)

# build tutorial (TODO)
# INCLUDE(UseLATEX.cmake)
