add_executable(Es4Opt Es4Opt.c)
target_link_libraries(Es4Opt PRIVATE cp)
set_target_properties(Es4Opt PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/programs/linear_systems/Es4Opt")

