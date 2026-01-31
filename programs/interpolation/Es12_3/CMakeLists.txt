add_executable(Es12_3 Es12_3.c)
target_link_libraries(Es12_3 PRIVATE cp)
set_target_properties(Es12_3 PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/programs/interpolation/Es12_3")

