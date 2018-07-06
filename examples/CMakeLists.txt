#
# Make the C++ example
#

add_executable(test_cplantbox main.cpp)
target_link_libraries(test_cplantbox CPlantBox)

# link the model parameter folder (on Windows we have to copy since symlinks are not supported)
if(${CMAKE_SYSTEM_NAME} STREQUAL "Windows")
    execute_process(COMMAND ${CMAKE_COMMAND} "-E" "copy" "${PROJECT_SOURCE_DIR}/modelparameter" "${CMAKE_CURRENT_BINARY_DIR}/modelparameter")
else()
    execute_process(COMMAND ${CMAKE_COMMAND} "-E" "create_symlink" "${PROJECT_SOURCE_DIR}/modelparameter" "${CMAKE_CURRENT_BINARY_DIR}/modelparameter")
endif()
