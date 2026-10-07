add_executable(true_anomaly true_anomaly.c)
target_link_libraries(true_anomaly PRIVATE cp)
set_target_properties(true_anomaly PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/programs/root_finding/es4_2_3_true_anomaly_calculation")

