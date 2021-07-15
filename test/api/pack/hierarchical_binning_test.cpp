#include <gtest/gtest.h>

#include <fstream>
#include <sstream>
#include <vector>

#include <chopper/pack/hierarchical_binning.hpp>

#include "../api_test.hpp"

TEST(hierarchical_binning_test, small_example)
{
    pack_config config;
    config.t_max = 4;

    std::stringstream output_buffer;
    std::stringstream header_buffer;
    pack_data data;
    data.output_buffer = &output_buffer;
    data.header_buffer = &header_buffer;
    data.filenames = {"seq0", "seq1", "seq2", "seq3", "seq4", "seq5", "seq6",  "seq7"};
    data.kmer_counts = {500, 1000, 500, 500, 500, 500, 500, 500};

    hierarchical_binning algo{data, config};
    EXPECT_EQ(algo.execute(), 3); // #HIGH_LEVEL_IBF max_bin_id:3

    std::string expected_file
    {
        "#MERGED_BIN_1 max_bin_id:0\n"
        "#MERGED_BIN_2 max_bin_id:0\n"
        "#MERGED_BIN_3 max_bin_id:0\n"
        "#FILES\tBIN_INDICES\tNUMBER_OF_BINS\tEST_MAX_TB_SIZES\n"
        "seq7\t0\t1\t500\n"
        "seq5\t1;0\t1;7\t1000;72\n"
        "seq6\t1;7\t1;57\t1000;9\n"
        "seq3\t2;0\t1;7\t1000;72\n"
        "seq4\t2;7\t1;57\t1000;9\n"
        "seq1\t3;0\t1;8\t2000;125\n"
        "seq0\t3;8\t1;4\t2000;125\n"
        "seq2\t3;12\t1;52\t2000;10\n"
    };

    EXPECT_EQ(header_buffer.str() + output_buffer.str(), expected_file);
}

TEST(hierarchical_binning_test, another_example)
{
    pack_config config;
    config.t_max = 5;

    std::stringstream output_buffer;
    std::stringstream header_buffer;
    pack_data data;
    data.output_buffer = &output_buffer;
    data.header_buffer = &header_buffer;
    data.filenames = {"seq0", "seq1", "seq2", "seq3", "seq4", "seq5", "seq6",  "seq7"};
    data.kmer_counts = {50, 1000, 1000, 50, 5, 10, 10, 5};

    hierarchical_binning algo{data, config};
    EXPECT_EQ(algo.execute(), 1); // #HIGH_LEVEL_IBF max_bin_id:1

    std::string expected_file
    {
        "#MERGED_BIN_0;0 max_bin_id:45\n"
        "#MERGED_BIN_0 max_bin_id:0\n"
        "#FILES\tBIN_INDICES\tNUMBER_OF_BINS\tEST_MAX_TB_SIZES\n"
        "seq6\t0;0;0\t1;1;45\t130;30;1\n"
        "seq5\t0;0;45\t1;1;6\t130;30;2\n"
        "seq7\t0;0;51\t1;1;3\t130;30;2\n"
        "seq4\t0;0;54\t1;1;10\t130;30;1\n"
        "seq0\t0;1\t1;2\t130;25\n"
        "seq3\t0;3\t1;2\t130;25\n"
        "seq2\t1\t2\t500\n"
        "seq1\t3\t2\t500\n"
    };

    EXPECT_EQ(header_buffer.str() + output_buffer.str(), expected_file);
}

TEST(hierarchical_binning_test, knuts_example)
{
    pack_config config;
    config.alpha = 1;
    config.t_max = 5;

    std::stringstream output_buffer;
    std::stringstream header_buffer;
    pack_data data;
    data.output_buffer = &output_buffer;
    data.header_buffer = &header_buffer;
    data.filenames = {"seq0", "seq1", "seq2", "seq3", "seq4"};
    data.kmer_counts = {60, 600, 1000, 800, 800};

    hierarchical_binning algo{data, config};
    EXPECT_EQ(algo.execute(), 1);

    std::string expected_file
    {
        "#MERGED_BIN_0 max_bin_id:0\n"
        "#FILES\tBIN_INDICES\tNUMBER_OF_BINS\tEST_MAX_TB_SIZES\n"
        "seq1\t0;0\t1;39\t660;16\n"
        "seq0\t0;39\t1;25\t660;3\n"
        "seq4\t1\t1\t800\n"
        "seq3\t2\t1\t800\n"
        "seq2\t3\t2\t500\n"
    };

    EXPECT_EQ(header_buffer.str() + output_buffer.str(), expected_file);
}
