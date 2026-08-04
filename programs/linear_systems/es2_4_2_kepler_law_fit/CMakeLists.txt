add_executable(keplerlawfit keplerlawfit.c)
target_link_libraries(keplerlawfit PRIVATE cp)
set_target_properties(keplerlawfit PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/programs/linear_systems/es2_4_2_kepler_law_fit")

