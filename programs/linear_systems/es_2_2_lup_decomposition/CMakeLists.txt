add_executable(Es6_1 Es6_1.c)
target_link_libraries(Es6_1 PRIVATE cp)
set_target_properties(Es6_1 PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/programs/linear_systems/Es6_1")

