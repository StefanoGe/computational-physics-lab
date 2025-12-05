set(EXERCISE_SOURCES Es2_2.c)
add_executable(Es2_2 ${EXERCISE_SOURCES})
target_link_libraries(Es2_2 PRIVATE cp_library m)
target_include_directories(Es2_2 PRIVATE ${CMAKE_SOURCE_DIR}/Libraries/header)
