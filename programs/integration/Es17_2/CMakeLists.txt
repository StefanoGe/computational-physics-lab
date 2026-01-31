add_executable(Es17_2 Es17_2.c)
target_link_libraries(Es17_2 PRIVATE cp)
set_target_properties(Es17_2 PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/programs/integration/Es17_2")

