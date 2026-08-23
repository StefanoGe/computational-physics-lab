add_executable(Es15_3 Es15_3.c)
target_link_libraries(Es15_3 PRIVATE cp)
set_target_properties(Es15_3 PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/programs/root_finding/Es15_3")

