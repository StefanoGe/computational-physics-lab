add_executable(Es8_2 Es8_2.c)
target_link_libraries(Es8_2 PRIVATE cp)
set_target_properties(Es8_2 PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/programs/linear_systems/Es8_2")

