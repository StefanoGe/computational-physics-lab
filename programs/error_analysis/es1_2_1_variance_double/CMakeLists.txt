add_executable(variancedouble variancedouble.c)
target_link_libraries(variancedouble PRIVATE cp)
set_target_properties(variancedouble PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/programs/error_analysis/es1_2_1_variance_double")

