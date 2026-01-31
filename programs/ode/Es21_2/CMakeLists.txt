add_executable(feuler_second_order feuler_second_order.c)
target_link_libraries(feuler_second_order PRIVATE cp)
set_target_properties(feuler_second_order PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/programs/ode/Es21_2")

