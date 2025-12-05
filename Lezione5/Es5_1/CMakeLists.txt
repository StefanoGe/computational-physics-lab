set(EXERCISE_SOURCES Es5_1.c)
add_executable(Es5_1 ${EXERCISE_SOURCES})
target_link_libraries(Es5_1 PRIVATE cp_library m)
target_include_directories(Es5_1 PRIVATE ${CMAKE_SOURCE_DIR}/Libraries/header)
