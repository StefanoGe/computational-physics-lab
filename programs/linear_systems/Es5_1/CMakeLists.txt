add_executable(Es5_1 Es5_1.c)
target_link_libraries(Es5_1 PRIVATE cp)
set_target_properties(Es5_1 PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/programs/linear_systems/Es5_1")

