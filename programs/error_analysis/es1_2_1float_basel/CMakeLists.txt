add_executable(baselfloat baselfloat.c)
target_link_libraries(baselfloat PRIVATE cp)
set_target_properties(baselfloat PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/programs/error_analysis/es1_2_1float_basel")

