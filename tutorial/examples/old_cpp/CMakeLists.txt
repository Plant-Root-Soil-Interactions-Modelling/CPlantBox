add_executable(test test.cpp)
set_target_properties(test PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/)
target_link_libraries(test CPlantBox)