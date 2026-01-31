add_executable(Es17_1 Es17_1.c)
target_link_libraries(Es17_1 PRIVATE cp)
set_target_properties(Es17_1 PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/programs/integration/Es17_1")

