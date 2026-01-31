add_executable(runge_kutta_test runge_kutta_test.c)
target_link_libraries(runge_kutta_test PRIVATE cp)
set_target_properties(runge_kutta_test PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/programs/ode/Es22_1")

