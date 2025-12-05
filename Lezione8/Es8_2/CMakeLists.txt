set(EXERCISE_SOURCES Es8_2.c)
add_executable(Es8_2 ${EXERCISE_SOURCES})
target_link_libraries(Es8_2 PRIVATE cp_library m)
target_include_directories(Es8_2 PRIVATE ${CMAKE_SOURCE_DIR}/Libraries/header)
