add_executable(Es12_2 Es12_2.c)
target_link_libraries(Es12_2 PRIVATE cp)
set_target_properties(Es12_2 PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/programs/interpolation/Es12_2")

