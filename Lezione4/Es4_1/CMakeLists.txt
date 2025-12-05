set(EXERCISE_SOURCES Es4_1.c)
add_executable(Es4_1 ${EXERCISE_SOURCES})
target_link_libraries(Es4_1 PRIVATE cp_library m)
target_include_directories(Es4_1 PRIVATE ${CMAKE_SOURCE_DIR}/Libraries/header)
