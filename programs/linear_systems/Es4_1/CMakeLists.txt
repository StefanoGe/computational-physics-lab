add_executable(Es4_1 Es4_1.c)
target_link_libraries(Es4_1 PRIVATE cp)
set_target_properties(Es4_1 PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/programs/linear_systems/Es4_1")

