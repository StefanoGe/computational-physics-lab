add_executable(Es4_1lapack Es4_1lapack.c)
target_link_libraries(Es4_1lapack PRIVATE cp)
set_target_properties(Es4_1lapack PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/programs/linear_systems/es2_1_1lapack_forw_backsubstitution")

