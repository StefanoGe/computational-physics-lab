add_executable(schroedinger_gaussian1 schroedinger_gaussian1.c)
target_link_libraries(schroedinger_gaussian1 PRIVATE cp)
set_target_properties(schroedinger_gaussian1 PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/programs/ode/Es23_1")

