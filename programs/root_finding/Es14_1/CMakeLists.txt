add_executable(Es14_1 Es14_1.c)
target_link_libraries(Es14_1 PRIVATE cp)
set_target_properties(Es14_1 PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/programs/root_finding/Es14_1")

