add_executable(poschl_teller poschl_teller.c)
target_link_libraries(poschl_teller PRIVATE cp)
set_target_properties(poschl_teller PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/programs/schroedinger/Es23_4")

