add_executable(simpson_rule simpson_rule.c)
target_link_libraries(simpson_rule PRIVATE cp)
set_target_properties(simpson_rule PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/programs/integration/es5_1_2_simpson_rule")

