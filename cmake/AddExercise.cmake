function(add_exercise ex_dir)
  # ex_dir like error_analysis/Ex1_1
  # You can change naming convention here later
  string(REPLACE "/" "_" target "${ex_dir}")   # e.g. error_analysis_Ex1_1

  add_executable(${target} ${CMAKE_CURRENT_SOURCE_DIR}/main.c)
  target_link_libraries(${target} PRIVATE cp)
endfunction()

