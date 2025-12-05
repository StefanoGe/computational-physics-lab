set(EXERCISE_SOURCES Es15_2.c)
add_executable(Es15_2 ${EXERCISE_SOURCES})
target_link_libraries(Es15_2 PRIVATE cp_library m)
target_include_directories(Es15_2 PRIVATE ${CMAKE_SOURCE_DIR}/Libraries/header)
