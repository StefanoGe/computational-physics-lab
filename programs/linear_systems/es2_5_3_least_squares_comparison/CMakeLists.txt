add_executable(Es9_2 Es9_2.c)
target_link_libraries(Es9_2 PRIVATE cp)
set_target_properties(Es9_2 PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/programs/linear_systems/Es9_2")

