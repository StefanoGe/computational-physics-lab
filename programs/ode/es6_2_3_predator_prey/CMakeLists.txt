add_executable(predator_prey predator_prey.c)
target_link_libraries(predator_prey PRIVATE cp)
set_target_properties(predator_prey PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/programs/ode/es6_2_3_predator_prey")

