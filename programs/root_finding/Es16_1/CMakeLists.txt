add_executable(Es16_1 Es16_1.c)
target_link_libraries(Es16_1 PRIVATE cp)
set_target_properties(Es16_1 PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/programs/root_finding/Es16_1")

