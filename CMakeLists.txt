cmake_minimum_required(VERSION 3.10)

project(ComputationalPhysicsLaboratory C)

# c standard
set(CMAKE_C_STANDARD 23)
set(CMAKE_C_STANDARD_REQUIRED ON)

# common flags
add_compile_options(-Wall -Wextra -Werror -pedantic-errors -Wconversion -Wno-sign-conversion -O0 -lm)

add_subdirectory(libs)
add_subdirectory(programs)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
