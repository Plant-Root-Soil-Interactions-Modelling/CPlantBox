cmake_minimum_required (VERSION 2.8.12)
project (CPlantBox)

# set compiler options
set(RELEASE_WARNING_OPTS "-Wall -Wunused -Wmissing-include-dirs -Wcast-align -Wno-sign-compare  ")
set(DEBUG_WARNING_OPTS "-Wall -Wunused -Wmissing-include-dirs -Wcast-align -Wno-sign-compare ")
set(RELEASE_OPTS "-std=c++11 -fno-strict-aliasing -fstrict-overflow -fno-finite-math-only -O3 -march=native -funroll-loops -g0 ")
set(CMAKE_CXX_FLAGS_RELEASE "${RELEASE_WARNING_OPTS} ${RELEASE_OPTS}")
set(CMAKE_CXX_FLAGS_DEBUG "-std=c++11 -O0 -ggdb -Wall")

# set default build type to release Release
set(CMAKE_BUILD_TYPE Release)

#compile options 
#to use photosynthesis:
#cmake . -DWITH_PHOTOSYNTHESIS=ON
#make
#to use photosynthesis + phloem:
#cmake . -DWITH_PHLOEM=ON
#make

option(WITH_PHLOEM "Option description" OFF)
option(WITH_PHOTOSYNTHESIS "Option description" OFF)

if (WITH_PHLOEM)
MESSAGE(STATUS "WITH_PHLOEM")
add_definitions(-DUSE_PHLOEM=1)
    add_definitions(-DUSE_PHOTOSYNTHESIS=1)
endif()

if (WITH_PHOTOSYNTHESIS)
MESSAGE(STATUS "WITH_PHOTOSYNTHESIS")
#need C++17 for Eigen3 
set(RELEASE_OPTS "-std=c++17 -fno-strict-aliasing -fstrict-overflow -fno-finite-math-only -O3 -march=native -funroll-loops -g0 ")
set(CMAKE_CXX_FLAGS_RELEASE "${RELEASE_WARNING_OPTS} ${RELEASE_OPTS}")
set(CMAKE_CXX_FLAGS_DEBUG "-std=c++17 -O0 -ggdb -Wall")

add_definitions(-DUSE_PHOTOSYNTHESIS=1)
find_package(Eigen3 REQUIRED CONFIG)

# checking property of the target
get_target_property(inc_dir Eigen3::Eigen INTERFACE_INCLUDE_DIRECTORIES)
message("Eigen include directory: ${inc_dir}")
include_directories(${inc_dir} )
endif()

# add source directory to the include path
include_directories(${PROJECT_SOURCE_DIR}/src)

# add subdirectories
add_subdirectory(src)
add_subdirectory(tutorial)