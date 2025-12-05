set(EXERCISE_SOURCES Es17_2.c)
add_executable(Es17_2 ${EXERCISE_SOURCES})
target_link_libraries(Es17_2 PRIVATE cp_library m)
target_include_directories(Es17_2 PRIVATE ${CMAKE_SOURCE_DIR}/Libraries/header)
