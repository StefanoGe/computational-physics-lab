add_executable(improper_integrals improper_integrals.c)
target_link_libraries(improper_integrals PRIVATE cp)
set_target_properties(improper_integrals PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/programs/integration/Es20_1")

