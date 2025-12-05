set(EXERCISE_SOURCES Es8_1.c)
add_executable(Es8_1 ${EXERCISE_SOURCES})
target_link_libraries(Es8_1 PRIVATE cp_library m)
target_include_directories(Es8_1 PRIVATE ${CMAKE_SOURCE_DIR}/Libraries/header)
