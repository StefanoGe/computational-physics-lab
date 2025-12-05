set(EXERCISE_SOURCES Es4Opt.c)
add_executable(Es4Opt ${EXERCISE_SOURCES})
target_link_libraries(Es4Opt PRIVATE cp_library m)
target_include_directories(Es4Opt PRIVATE ${CMAKE_SOURCE_DIR}/Libraries/header)
