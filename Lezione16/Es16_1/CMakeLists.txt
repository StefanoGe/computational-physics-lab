set(EXERCISE_SOURCES Es16_1.c)
add_executable(Es16_1 ${EXERCISE_SOURCES})
target_link_libraries(Es16_1 PRIVATE cp_library m)
target_include_directories(Es16_1 PRIVATE ${CMAKE_SOURCE_DIR}/Libraries/header)
