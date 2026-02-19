cmake_minimum_required(VERSION 3.25)
project(CPlantBox LANGUAGES CXX)


include(CMakePackageConfigHelpers) # for later package export

# use static CRT on MSVC to avoid needing the VC++ redistributable
#if (MSVC)
#  set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
#endif()

#default for single-config generators
if (NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
endif()

# subproject
add_subdirectory(src)