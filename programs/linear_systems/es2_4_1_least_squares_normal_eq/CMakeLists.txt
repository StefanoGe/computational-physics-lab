add_executable(Es8_1 Es8_1.c)
target_link_libraries(Es8_1 PRIVATE cp)
set_target_properties(Es8_1 PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/programs/linear_systems/es2_4_1_least_squares_normal_eq")

