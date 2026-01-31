add_executable(variancefloat variancefloat.c)
target_link_libraries(variancefloat PRIVATE cp)
set_target_properties(variancefloat PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/programs/error_analysis/Es3_1float")

