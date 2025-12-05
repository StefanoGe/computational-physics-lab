set(EXERCISE_SOURCES Es12_1.c)
add_executable(Es12_1 ${EXERCISE_SOURCES})
target_link_libraries(Es12_1 PRIVATE cp_library m)
target_include_directories(Es12_1 PRIVATE ${CMAKE_SOURCE_DIR}/Libraries/header)
