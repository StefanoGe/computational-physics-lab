set(EXERCISE_SOURCES Es12_3.c)
add_executable(Es12_3 ${EXERCISE_SOURCES})
target_link_libraries(Es12_3 PRIVATE cp_library m)
target_include_directories(Es12_3 PRIVATE ${CMAKE_SOURCE_DIR}/Libraries/header)
