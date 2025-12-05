set(EXERCISE_SOURCES Es3_1float.c)
add_executable(Es3_1float ${EXERCISE_SOURCES})
target_link_libraries(Es3_1float PRIVATE cp_library m)
target_include_directories(Es3_1float PRIVATE ${CMAKE_SOURCE_DIR}/Libraries/header)
