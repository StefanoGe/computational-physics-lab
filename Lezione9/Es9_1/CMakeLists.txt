set(EXERCISE_SOURCES Es9_1.c)
add_executable(Es9_1 ${EXERCISE_SOURCES})
target_link_libraries(Es9_1 PRIVATE cp_library m)
target_include_directories(Es9_1 PRIVATE ${CMAKE_SOURCE_DIR}/Libraries/header)
