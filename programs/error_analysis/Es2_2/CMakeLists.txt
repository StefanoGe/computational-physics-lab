add_executable(sequence sequence.c)
target_link_libraries(sequence PRIVATE cp)
set_target_properties(sequence PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/programs/error_analysis/Es2_2")

