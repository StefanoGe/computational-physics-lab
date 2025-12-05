set(EXERCISE_SOURCES Es7_1.c)
add_executable(Es7_1 ${EXERCISE_SOURCES})
target_link_libraries(Es7_1 PRIVATE cp_library m)
target_include_directories(Es7_1 PRIVATE ${CMAKE_SOURCE_DIR}/Libraries/header)
