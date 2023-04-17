# don't build per default

add_executable(test test_plant.cpp)
set_target_properties(test PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/)
include_directories(${PROJECT_BINARY_DIR}/src/external/tinyxml2)
target_link_libraries(test CPlantBox)
