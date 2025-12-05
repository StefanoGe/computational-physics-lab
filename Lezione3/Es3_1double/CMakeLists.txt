set(EXERCISE_SOURCES Es3_1double.c)
add_executable(Es3_1double ${EXERCISE_SOURCES})
target_link_libraries(Es3_1double PRIVATE cp_library m)
target_include_directories(Es3_1double PRIVATE ${CMAKE_SOURCE_DIR}/Libraries/header)
