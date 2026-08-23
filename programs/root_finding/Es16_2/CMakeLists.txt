add_executable(Es16_2 Es16_2.c)
target_link_libraries(Es16_2 PRIVATE cp)
set_target_properties(Es16_2 PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/programs/root_finding/Es16_2")

