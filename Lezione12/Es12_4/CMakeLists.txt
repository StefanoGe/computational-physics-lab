set(EXERCISE_SOURCES Es12_4.c)
add_executable(Es12_4 ${EXERCISE_SOURCES})
target_link_libraries(Es12_4 PRIVATE cp_library m)
target_include_directories(Es12_4 PRIVATE ${CMAKE_SOURCE_DIR}/Libraries/header)
