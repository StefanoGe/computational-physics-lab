add_executable(feuler_first_order feuler_first_order.c)
target_link_libraries(feuler_first_order PRIVATE cp)
set_target_properties(feuler_first_order PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/programs/ode/Es21_1")

