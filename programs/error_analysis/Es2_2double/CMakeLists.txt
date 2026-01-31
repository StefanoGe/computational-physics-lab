add_executable(baseldouble baseldouble.c)
target_link_libraries(baseldouble PRIVATE cp)
set_target_properties(baseldouble PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/programs/error_analysis/Es2_2double")

