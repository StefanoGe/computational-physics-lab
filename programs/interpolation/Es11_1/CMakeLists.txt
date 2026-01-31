add_executable(Es11_1 Es11_1.c)
target_link_libraries(Es11_1 PRIVATE cp)
set_target_properties(Es11_1 PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/programs/interpolation/Es11_1")

