add_executable(Es10_1 Es10_1.c)
target_link_libraries(Es10_1 PRIVATE cp)
set_target_properties(Es10_1 PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/programs/interpolation/Es10_1")

