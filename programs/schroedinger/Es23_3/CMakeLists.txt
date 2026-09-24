add_executable(finite_barrier finite_barrier.c)
target_link_libraries(finite_barrier PRIVATE cp)
set_target_properties(finite_barrier PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/programs/schroedinger/Es23_3")

