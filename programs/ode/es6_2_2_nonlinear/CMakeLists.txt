add_executable(rk_nonlinear rk_nonlinear.c)
target_link_libraries(rk_nonlinear PRIVATE cp)
set_target_properties(rk_nonlinear PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/programs/ode/es6_2_2_nonlinear")

