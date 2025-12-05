cmake_minimum_required(VERSION 3.10)

project(ComputationalPhysicsLaboratory C)

# c standard
set(CMAKE_C_STANDARD 23)
set(CMAKE_C_STANDARD_REQUIRED ON)

# common flags
add_compile_options(-Wall -Wextra -Werror -pedantic-errors -Wconversion -Wno-sign-conversion -O0 -lm)

set(LAB_BIN_DIR ${CMAKE_SOURCE_DIR}/bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${LAB_BIN_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${LAB_BIN_DIR})

# include headers
include_directories(${CMAKE_SOURCE_DIR}/Libraries/header)

# make library
file(GLOB LIB_SOURCES "${CMAKE_SOURCE_DIR}/Libraries/source/*.c")
add_library(cp_library STATIC ${LIB_SOURCES})
target_include_directories(cp_library PUBLIC ${CMAKE_SOURCE_DIR}/Libraries/header)

add_subdirectory(Lezione18/Es18_1)
add_subdirectory(LezioneTest/es_test/)
add_subdirectory(LezioneTest/test/)
