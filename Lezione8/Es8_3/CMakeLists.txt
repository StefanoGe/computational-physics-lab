set(EXERCISE_SOURCES Es8_3.c)
add_executable(Es8_3 ${EXERCISE_SOURCES})
target_link_libraries(Es8_3 PRIVATE cp_library m)
target_include_directories(Es8_3 PRIVATE ${CMAKE_SOURCE_DIR}/Libraries/header)
