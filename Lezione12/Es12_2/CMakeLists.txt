set(EXERCISE_SOURCES Es12_2.c)
add_executable(Es12_2 ${EXERCISE_SOURCES})
target_link_libraries(Es12_2 PRIVATE cp_library m)
target_include_directories(Es12_2 PRIVATE ${CMAKE_SOURCE_DIR}/Libraries/header)
