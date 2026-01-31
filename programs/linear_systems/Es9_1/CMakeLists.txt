add_executable(Es9_1 Es9_1.c)
target_link_libraries(Es9_1 PRIVATE cp)
set_target_properties(Es9_1 PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/programs/linear_systems/Es9_1")

