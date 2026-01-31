file(GLOB EX_DIRS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *)

foreach(ex ${EX_DIRS})
    if(IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${ex}
       AND EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${child}/CMakeLists.txt)
        add_subdirectory(${ex})
    endif()
endforeach()

