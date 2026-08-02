add_executable(qrtest qrtest.c)
target_link_libraries(qrtest PRIVATE cp)
set_target_properties(qrtest PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/programs/linear_systems/es2_5_1_qr_test")

