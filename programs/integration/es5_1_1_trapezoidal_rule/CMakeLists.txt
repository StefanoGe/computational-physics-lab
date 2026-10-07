add_executable(trapezoidal_rule trapezoidal_rule.c)
target_link_libraries(trapezoidal_rule PRIVATE cp)
set_target_properties(trapezoidal_rule PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/programs/integration/es5_1_1_trapezoidal_rule")

