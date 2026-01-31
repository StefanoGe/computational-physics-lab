file(GLOB CP_SOURCES src/*.c)

add_library(cp STATIC ${CP_SOURCES})

target_include_directories(cp PUBLIC include)
target_link_libraries(cp PUBLIC m)

