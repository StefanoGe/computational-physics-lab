set(EXERCISE_SOURCES Es17_1.c)
add_executable(Es17_1 ${EXERCISE_SOURCES})
target_link_libraries(Es17_1 PRIVATE cp_library m)
target_include_directories(Es17_1 PRIVATE ${CMAKE_SOURCE_DIR}/Libraries/header)
