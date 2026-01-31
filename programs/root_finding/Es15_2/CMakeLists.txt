add_executable(Es15_2 Es15_2.c)
target_link_libraries(Es15_2 PRIVATE cp)
set_target_properties(Es15_2 PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/programs/root_finding/Es15_2")

