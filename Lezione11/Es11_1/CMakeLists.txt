set(EXERCISE_SOURCES Es11_1.c)
add_executable(Es11_1 ${EXERCISE_SOURCES})
target_link_libraries(Es11_1 PRIVATE cp_library m)
target_include_directories(Es11_1 PRIVATE ${CMAKE_SOURCE_DIR}/Libraries/header)
