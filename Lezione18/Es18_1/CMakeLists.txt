set(EXERCISE_SOURCES gauss_legendre_weights.c)
add_executable(gauss_legendre_weights ${EXERCISE_SOURCES})
target_link_libraries(gauss_legendre_weights PRIVATE cp_library m)
target_include_directories(gauss_legendre_weights PRIVATE ${CMAKE_SOURCE_DIR}/Libraries/header)
