add_executable(tise tise.c)
target_link_libraries(tise PRIVATE cp)
set_target_properties(tise PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/programs/schroedinger/es7_2_1_tise")

