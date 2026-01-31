add_executable(lotka_volterra lotka_volterra.c)
target_link_libraries(lotka_volterra PRIVATE cp)
set_target_properties(lotka_volterra PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/programs/ode/Es22_3")

