add_executable(Es3_1double Es3_1double.c)
target_link_libraries(Es3_1double PRIVATE cp)
set_target_properties(Es3_1double PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/programs/error_analysis/Es3_1double")

