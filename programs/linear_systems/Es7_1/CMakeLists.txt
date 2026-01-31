add_executable(Es7_1 Es7_1.c)
target_link_libraries(Es7_1 PRIVATE cp)
set_target_properties(Es7_1 PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/programs/linear_systems/Es7_1")

