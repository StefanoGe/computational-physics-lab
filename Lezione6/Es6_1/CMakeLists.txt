set(EXERCISE_SOURCES Es6_1.c)
add_executable(Es6_1 ${EXERCISE_SOURCES})
target_link_libraries(Es6_1 PRIVATE cp_library m)
target_include_directories(Es6_1 PRIVATE ${CMAKE_SOURCE_DIR}/Libraries/header)
