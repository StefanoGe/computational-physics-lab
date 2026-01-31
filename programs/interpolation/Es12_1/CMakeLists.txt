add_executable(Es12_1 Es12_1.c)
target_link_libraries(Es12_1 PRIVATE cp)
set_target_properties(Es12_1 PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/programs/interpolation/Es12_1")

