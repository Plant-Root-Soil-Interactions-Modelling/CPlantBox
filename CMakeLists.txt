cmake_minimum_required (VERSION 3.5)
project (CPlantBox LANGUAGES CXX)

# set compiler options
set(RELEASE_WARNING_OPTS "-Wall -Wunused -Wmissing-include-dirs -Wcast-align -Wno-sign-compare  ")
set(DEBUG_WARNING_OPTS "-Wall -Wunused -Wmissing-include-dirs -Wcast-align -Wno-sign-compare ")
# Portable release opts (no -march=native for wheel portability)
set(RELEASE_OPTS "-std=c++17 -fno-strict-aliasing -fstrict-overflow -fno-finite-math-only -O3 -funroll-loops -g0 ")
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

# CTest smoke hook to run curated pytest subset from source root
include(CTest)
if(BUILD_TESTING)
  add_test(
    NAME smoke_pytest
    COMMAND python3 -m pytest -q 
            test/test_organ.py 
            test/test_organparameter.py 
            test/test_rootparameter.py 
            test/test_leafparameter.py 
            test/test_leaf.py 
            test/test_hydraulic_parameters.py
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  )
  # Deterministic defaults and module path
  set_tests_properties(smoke_pytest PROPERTIES 
    ENVIRONMENT "OMP_NUM_THREADS=1;PYTHONPATH=${CMAKE_SOURCE_DIR};${CMAKE_SOURCE_DIR}/src"
  )
endif()
