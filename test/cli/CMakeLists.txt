cmake_minimum_required (VERSION 3.18)

add_cli_test (cli_chopper_basic_test.cpp)

add_cli_test (cli_chopper_pipeline_test.cpp)
target_use_datasources (cli_chopper_pipeline_test FILES seq1.fa)
target_use_datasources (cli_chopper_pipeline_test FILES small.fa)
target_use_datasources (cli_chopper_pipeline_test FILES small2.fa)
target_use_datasources (cli_chopper_pipeline_test FILES small.split)
