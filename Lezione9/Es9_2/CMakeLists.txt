set(EXERCISE_SOURCES Es9_2.c)
add_executable(Es9_2 ${EXERCISE_SOURCES})
target_link_libraries(Es9_2 PRIVATE cp_library m)
target_include_directories(Es9_2 PRIVATE ${CMAKE_SOURCE_DIR}/Libraries/header)
