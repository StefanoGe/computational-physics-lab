add_executable(eigenvalues eigenvalues.c)
target_link_libraries(eigenvalues PRIVATE cp)
set_target_properties(eigenvalues PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/programs/linear_systems/es2_5_2_eigenvalues")

