add_executable(Es12_4 Es12_4.c)
target_link_libraries(Es12_4 PRIVATE cp)
set_target_properties(Es12_4 PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/programs/interpolation/Es12_4")

