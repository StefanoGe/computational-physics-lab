add_executable(vanpolcosine vanpolcosine.c)
target_link_libraries(vanpolcosine PRIVATE cp)
set_target_properties(vanpolcosine PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/programs/interpolation/es3_1_1_polynomial_int_vandermonde")

