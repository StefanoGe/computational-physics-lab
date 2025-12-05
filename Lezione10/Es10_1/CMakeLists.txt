set(EXERCISE_SOURCES Es10_1.c)
add_executable(Es10_1 ${EXERCISE_SOURCES})
target_link_libraries(Es10_1 PRIVATE cp_library m)
target_include_directories(Es10_1 PRIVATE ${CMAKE_SOURCE_DIR}/Libraries/header)
