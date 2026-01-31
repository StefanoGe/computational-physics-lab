add_executable(Es15_1 Es15_1.c)
target_link_libraries(Es15_1 PRIVATE cp)
set_target_properties(Es15_1 PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/programs/root_finding/Es15_1")

