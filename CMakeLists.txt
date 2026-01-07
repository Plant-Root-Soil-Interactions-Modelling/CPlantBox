cmake_minimum_required (VERSION 3.5)
project (CPlantBox LANGUAGES CXX)

# set compiler options
set(RELEASE_WARNING_OPTS "-Wall -Wunused -Wmissing-include-dirs -Wcast-align -Wno-sign-compare  ")
set(DEBUG_WARNING_OPTS "-Wall -Wunused -Wmissing-include-dirs -Wcast-align -Wno-sign-compare ")
set(RELEASE_OPTS "-std=c++17 -fno-strict-aliasing -fstrict-overflow -fno-finite-math-only -O3 -march=native -funroll-loops -g0 ")
set(CMAKE_CXX_FLAGS_RELEASE "${RELEASE_WARNING_OPTS} ${RELEASE_OPTS}")
set(CMAKE_CXX_FLAGS_DEBUG "-std=c++17 -O0 -ggdb -Wall")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# set default build type to release Release, Debug, RelWithDebInfo 
set(CMAKE_BUILD_TYPE RelWithDebInfo)

# add subdirectories
add_subdirectory(src)
#add_subdirectory(tutorial)

# build tutorial (TODO)
# INCLUDE(UseLATEX.cmake)
