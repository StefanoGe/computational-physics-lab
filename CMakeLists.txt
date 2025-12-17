cmake_minimum_required(VERSION 3.10)

project(ComputationalPhysicsLaboratory C)

# c standard
set(CMAKE_C_STANDARD 23)
set(CMAKE_C_STANDARD_REQUIRED ON)

# common flags
add_compile_options(-Wall -Wextra -Werror -pedantic-errors -Wconversion -Wno-sign-conversion -O0 -lm)

#set(LAB_BIN_DIR ${CMAKE_SOURCE_DIR}/bin)
#set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${LAB_BIN_DIR})

# include headers
include_directories(${CMAKE_SOURCE_DIR}/Libraries/header)

# make library
file(GLOB LIB_SOURCES "${CMAKE_SOURCE_DIR}/Libraries/source/*.c")
add_library(cp_library STATIC ${LIB_SOURCES})
target_include_directories(cp_library PUBLIC ${CMAKE_SOURCE_DIR}/Libraries/header)

add_subdirectory(LezioneTest/es_test/)
add_subdirectory(LezioneTest/test/)
add_subdirectory(Lezione10/Es10_1/)
add_subdirectory(Lezione11/Es11_1/)
add_subdirectory(Lezione12/Es12_1/)
add_subdirectory(Lezione12/Es12_2/)
add_subdirectory(Lezione12/Es12_3/)
add_subdirectory(Lezione12/Es12_4/)
add_subdirectory(Lezione15/Es15_1/)
add_subdirectory(Lezione15/Es15_2/)
add_subdirectory(Lezione16/Es16_1/)
add_subdirectory(Lezione17/Es17_1/)
add_subdirectory(Lezione17/Es17_2/)
add_subdirectory(Lezione18/Es18_1/)
add_subdirectory(Lezione2/Es2_1/)
add_subdirectory(Lezione2/Es2_2/)
add_subdirectory(Lezione3/Es3_1double/)
add_subdirectory(Lezione3/Es3_1float/)
add_subdirectory(Lezione4/Es4_1/)
add_subdirectory(Lezione4/Es4Opt/)
add_subdirectory(Lezione5/Es5_1/)
add_subdirectory(Lezione6/Es6_1/)
add_subdirectory(Lezione7/Es7_1/)
add_subdirectory(Lezione8/Es8_1/)
add_subdirectory(Lezione8/Es8_2/)
add_subdirectory(Lezione8/Es8_3/)
add_subdirectory(Lezione9/Es9_1/)
add_subdirectory(Lezione9/Es9_2/)
add_subdirectory(Lezione19/Es19_1/)
add_subdirectory(Lezione14/Es14_1/)
add_subdirectory(Lezione1/Es1_1/)
add_subdirectory(Lezione20/Es20_1/)
add_subdirectory(Lezione21/Es21_1/)
add_subdirectory(Lezione21/Es21_2/)
add_subdirectory(Lezione22/Es22_1/)
add_subdirectory(Lezione22/rk_nonlinear/)
add_subdirectory(Lezione23/Es23_1/)
