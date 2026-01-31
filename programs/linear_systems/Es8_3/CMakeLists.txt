add_executable(Es8_3 Es8_3.c)
target_link_libraries(Es8_3 PRIVATE cp)
set_target_properties(Es8_3 PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/programs/linear_systems/Es8_3")

