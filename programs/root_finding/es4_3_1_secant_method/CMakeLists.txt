add_executable(secant_method secant_method.c)
target_link_libraries(secant_method PRIVATE cp)
set_target_properties(secant_method PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/programs/root_finding/es4_3_1_secant_method")

