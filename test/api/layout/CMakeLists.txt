cmake_minimum_required (VERSION 3.18)

add_api_test (ibf_query_cost_test.cpp)
add_api_test (execute_layout_test.cpp)

add_api_test (execute_with_estimation_test.cpp)
target_use_datasources (execute_with_estimation_test FILES seq1.fa)
target_use_datasources (execute_with_estimation_test FILES seq2.fa)
target_use_datasources (execute_with_estimation_test FILES seq3.fa)
target_use_datasources (execute_with_estimation_test FILES small.fa)

add_api_test (fp_correction_test.cpp)
add_api_test (hibf_statistics_test.cpp)
add_api_test (layout_test.cpp)
