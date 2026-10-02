// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#include <gtest/gtest.h> // for Test, TestInfo, EXPECT_EQ, TEST

#include <algorithm> // for sort
#include <cstddef>   // for size_t
#include <cstdint>   // for uint64_t
#include <fstream>   // for ifstream
#include <numeric>   // for iota
#include <stdexcept> // for invalid_argument
#include <string>    // for string
#include <vector>    // for vector

#include <chopper/configuration.hpp>
#include <chopper/layout/execute.hpp>
#include <chopper/layout/input.hpp>
#include <chopper/layout/phibf/partition_user_bins.hpp>

#include <hibf/sketch/compute_sketches.hpp>
#include <hibf/sketch/hyperloglog.hpp>
#include <hibf/sketch/minhashes.hpp>

#include "../api_test.hpp"

namespace
{

// splitmix64 finaliser: turns consecutive integers into well-distributed hashes. HyperLogLog and the MinHash buckets
// (hash & 15, each needs 40 values) both assume uniformly distributed hashes, which plain consecutive integers are not.
uint64_t scramble(uint64_t x)
{
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

// 240 user bins in 24 families of 10 user bins. User bins of the same family share their first
// min(kmer_counts) hashes. Sizes range from 1'000 to 24'500 hashes.
struct phibf_data
{
    static constexpr size_t number_of_user_bins{240};

    std::vector<size_t> kmer_counts{};
    std::vector<size_t> content_ids{};
    chopper::configuration config{};
    std::vector<seqan::hibf::sketch::hyperloglog> sketches{};
    std::vector<seqan::hibf::sketch::minhashes> minHash_sketches{};
    std::vector<size_t> cardinalities{};

    phibf_data(size_t const number_of_partitions, int const partitioning_approach)
    {
        for (size_t ub = 0; ub < number_of_user_bins; ++ub)
        {
            content_ids.push_back(ub / 10);
            kmer_counts.push_back(1'000 + 500 * ((ub * 7) % 48));
        }

        config.number_of_partitions = number_of_partitions;
        config.partitioning_approach = partitioning_approach;
        config.hibf_config.tmax = 64;
        config.hibf_config.number_of_user_bins = number_of_user_bins;
        config.hibf_config.input_fn = [this](size_t const ub, seqan::hibf::insert_iterator it)
        {
            uint64_t const offset = static_cast<uint64_t>(content_ids[ub]) << 32;
            for (uint64_t i = 0; i < kmer_counts[ub]; ++i)
                it = scramble(offset + i);
        };
        config.hibf_config.validate_and_set_defaults();

        seqan::hibf::sketch::compute_sketches(config.hibf_config, sketches, minHash_sketches);

        for (auto const & sketch : sketches)
            cardinalities.push_back(sketch.estimate());
    }
};

// The data of the PHIBF stress test: `number_of_user_bins` user bins in families of 5 that mostly overlap.
// dist: "uniform" (3'000 to 23'000 hashes), "equal" (5'000 hashes), or "skewed" (two user bins with 400'000
// hashes, the rest 3'000 to 3'600). Failures found by the stress test are added as regression tests below.
struct stress_data
{
    std::vector<size_t> kmer_counts{};
    chopper::configuration config{};
    std::vector<seqan::hibf::sketch::hyperloglog> sketches{};
    std::vector<seqan::hibf::sketch::minhashes> minHash_sketches{};
    std::vector<size_t> cardinalities{};

    stress_data(int const partitioning_approach,
                size_t const number_of_user_bins,
                size_t const number_of_partitions,
                std::string const & dist)
    {
        for (size_t i = 0; i < number_of_user_bins; ++i)
        {
            if (dist == "equal")
                kmer_counts.push_back(5000);
            else if (dist == "skewed")
                kmer_counts.push_back((i < 2) ? 400'000 : 3000 + 100 * (i % 7));
            else
                kmer_counts.push_back(3000 + (scramble(i) % 20000));
        }

        config.number_of_partitions = number_of_partitions;
        config.partitioning_approach = partitioning_approach;
        config.hibf_config.number_of_user_bins = number_of_user_bins;
        config.hibf_config.tmax = 64;
        config.hibf_config.input_fn = [this](size_t const ub, seqan::hibf::insert_iterator it)
        {
            uint64_t const offset = static_cast<uint64_t>(ub / 5) << 32;
            for (uint64_t i = 0; i < kmer_counts[ub]; ++i)
                it = scramble(offset + i + (ub % 5) * 97);
        };
        config.hibf_config.validate_and_set_defaults();

        seqan::hibf::sketch::compute_sketches(config.hibf_config, sketches, minHash_sketches);

        for (auto const & sketch : sketches)
            cardinalities.push_back(sketch.estimate());
    }

    std::vector<std::vector<size_t>> partition() const
    {
        std::vector<std::vector<size_t>> partitions(config.number_of_partitions);
        chopper::layout::phibf::partition_user_bins(config, cardinalities, sketches, minHash_sketches, partitions);
        return partitions;
    }
};

// Expects that every user bin in [0, number_of_user_bins) is assigned to exactly one partition.
void expect_each_user_bin_assigned_once(std::vector<std::vector<size_t>> const & partitions,
                                        size_t const number_of_user_bins)
{
    std::vector<size_t> assigned_user_bins{};
    for (auto const & partition : partitions)
        assigned_user_bins.insert(assigned_user_bins.end(), partition.begin(), partition.end());
    std::ranges::sort(assigned_user_bins);

    std::vector<size_t> expected_user_bins(number_of_user_bins);
    std::iota(expected_user_bins.begin(), expected_user_bins.end(), 0u);
    EXPECT_EQ(assigned_user_bins, expected_user_bins);
}

} // namespace

class phibf_partition_test : public ::testing::TestWithParam<int>
{};

TEST_P(phibf_partition_test, each_user_bin_is_assigned_to_exactly_one_partition)
{
    size_t const number_of_partitions{4};
    phibf_data data{number_of_partitions, GetParam()};

    std::vector<std::vector<size_t>> partitions(number_of_partitions);
    chopper::layout::phibf::partition_user_bins(data.config,
                                                data.cardinalities,
                                                data.sketches,
                                                data.minHash_sketches,
                                                partitions);

    ASSERT_EQ(partitions.size(), number_of_partitions);
    expect_each_user_bin_assigned_once(partitions, phibf_data::number_of_user_bins);
}

INSTANTIATE_TEST_SUITE_P(all_approaches,
                         phibf_partition_test,
                         ::testing::Values(chopper::layout::phibf::partitioning_scheme::blocked,
                                           chopper::layout::phibf::partitioning_scheme::sorted,
                                           chopper::layout::phibf::partitioning_scheme::folded,
                                           chopper::layout::phibf::partitioning_scheme::weighted_fold,
                                           chopper::layout::phibf::partitioning_scheme::similarity,
                                           chopper::layout::phibf::partitioning_scheme::lsh,
                                           chopper::layout::phibf::partitioning_scheme::lsh_sim),
                         [](::testing::TestParamInfo<int> const & info)
                         {
                             return std::to_string(info.param);
                         });

TEST(phibf_execute_test, writes_one_layout_per_partition)
{
    seqan3::test::tmp_directory tmp_dir{};
    std::filesystem::path const layout_file{tmp_dir.path() / "phibf.layout"};

    size_t const number_of_partitions{3};
    phibf_data data{number_of_partitions, chopper::layout::phibf::partitioning_scheme::lsh_sim};
    data.config.output_filename = layout_file;

    std::vector<std::vector<std::string>> filenames{};
    for (size_t ub = 0; ub < phibf_data::number_of_user_bins; ++ub)
        filenames.push_back({"ub" + std::to_string(ub) + ".fa"});

    chopper::layout::execute(data.config, filenames, data.sketches, data.minHash_sketches);

    std::ifstream layout_stream{layout_file};
    auto [read_filenames, read_config, layouts] = chopper::layout::read_layouts_file(layout_stream);

    EXPECT_EQ(read_filenames, filenames);
    EXPECT_EQ(read_config.number_of_partitions, number_of_partitions);
    EXPECT_EQ(read_config.partitioning_approach, chopper::layout::phibf::partitioning_scheme::lsh_sim);
    ASSERT_EQ(layouts.size(), number_of_partitions);

    // Each layout holds the user bins of one partition, with global user bin indices.
    std::vector<std::vector<size_t>> partitions{};
    for (auto const & layout : layouts)
    {
        EXPECT_FALSE(layout.user_bins.empty());
        partitions.emplace_back();
        for (auto const & user_bin : layout.user_bins)
            partitions.back().push_back(user_bin.idx);
    }
    expect_each_user_bin_assigned_once(partitions, phibf_data::number_of_user_bins);
}

TEST(phibf_partition_test, unknown_partitioning_approach)
{
    size_t const number_of_partitions{4};

    for (int const approach : {-1, 7})
    {
        phibf_data data{number_of_partitions, approach};

        std::vector<std::vector<size_t>> partitions(number_of_partitions);
        EXPECT_THROW(chopper::layout::phibf::partition_user_bins(data.config,
                                                                 data.cardinalities,
                                                                 data.sketches,
                                                                 data.minHash_sketches,
                                                                 partitions),
                     std::invalid_argument);
    }
}

// post_process_clusters asserted that the clusters after the first number_of_partitions ones are sorted by the
// cardinality of their id() instead of their largest user bin.
TEST(phibf_regression_test, lsh_post_process_clusters_sanity_check)
{
    for (auto const & [n, np] : {std::pair<size_t, size_t>{100, 3}, {257, 2}})
    {
        stress_data const data{chopper::layout::phibf::partitioning_scheme::lsh, n, np, "equal"};
        expect_each_user_bin_assigned_once(data.partition(), n);
    }
}

// With fewer user bins than partitions, lsh and lsh_sim read and wrote out of bounds (std::partial_sort) and the other
// approaches left partitions empty.
TEST(phibf_regression_test, fewer_user_bins_than_partitions)
{
    for (int approach = 0; approach <= chopper::layout::phibf::partitioning_scheme::lsh_sim; ++approach)
    {
        stress_data const data{approach, 1, 2, "uniform"};
        EXPECT_THROW(data.partition(), std::invalid_argument) << "approach " << approach;
    }
}

// Some approaches leave partitions empty, e.g. blocked if there are fewer blocks than partitions, or sorted with a few
// very large user bins. An empty partition crashed compute_layout.
TEST(phibf_regression_test, empty_partition)
{
    using chopper::layout::phibf::partitioning_scheme;

    for (auto const & [approach, n, np, dist] : {std::tuple<int, size_t, size_t, std::string>{partitioning_scheme::blocked, 100, 16, "uniform"},
                                                 {partitioning_scheme::sorted, 64, 8, "skewed"}})
    {
        stress_data const data{approach, n, np, dist};
        EXPECT_THROW(data.partition(), std::runtime_error) << "approach " << approach;
    }
}

// sorted and folded moved on to the next partition whenever a partition reached its target cardinality. With
// zero-cardinality user bins at the end, the next one was written past the last partition.
TEST(phibf_regression_test, zero_cardinality_user_bins)
{
    using chopper::layout::phibf::partitioning_scheme;

    for (int const approach : {partitioning_scheme::sorted, partitioning_scheme::folded})
    {
        stress_data data{approach, 5, 2, "equal"};
        data.cardinalities = {10, 10, 10, 10, 0};

        auto const partitions = data.partition();
        expect_each_user_bin_assigned_once(partitions, 5);
        // The zero-cardinality user bin goes to the last partition (sorted) or the last folded part (folded).
        EXPECT_EQ(partitions[approach == partitioning_scheme::sorted ? 1 : 0].back(), 4u) << "approach " << approach;
    }
}

// similarity initialises each partition with one block of user bins, but the block size could result in fewer
// blocks than partitions, e.g. 33 user bins and 8 partitions: block size 5, 7 blocks.
TEST(phibf_regression_test, similarity_fewer_blocks_than_partitions)
{
    for (auto const & [n, np] : {std::pair<size_t, size_t>{33, 8}, {33, 16}, {100, 16}})
    {
        stress_data const data{chopper::layout::phibf::partitioning_scheme::similarity, n, np, "uniform"};
        expect_each_user_bin_assigned_once(data.partition(), n);
    }
}

// lsh and lsh_sim initialise each partition with one cluster. If there are fewer clusters than partitions, lsh_sim read
// the first cluster of an empty MultiCluster (assertion in Debug, out of bounds in Release). Now, the remaining
// partitions stay empty and partition_user_bins reports that.
TEST(phibf_regression_test, fewer_clusters_than_partitions)
{
    using chopper::layout::phibf::partitioning_scheme;

    for (int const approach : {partitioning_scheme::lsh, partitioning_scheme::lsh_sim})
    {
        for (size_t const n : {9u, 10u})
        {
            stress_data const data{approach, n, 8, "equal"};
            EXPECT_THROW(data.partition(), std::runtime_error) << "approach " << approach << " n " << n;
        }
    }
}

// find_best_partition only assigns a user bin to a partition below the target cardinality per partition that has room
// for it (up to 1.2 times the target). lsh and lsh_sim asserted that a single user bin always fits and otherwise left
// it unassigned in Release, which failed the final sanity check with a generic error. similarity ignored the result.
TEST(phibf_regression_test, user_bin_does_not_fit)
{
    using chopper::layout::phibf::partitioning_scheme;

    // similarity: 4 equally sized user bins (cardinality c) in 3 partitions. Each partition is initialised with one
    // user bin. The fourth would need room for 2c, but a partition may hold at most 1.2 * ceil(4c / 3) = 1.6c.
    for (auto const & [approach, n, np, dist] : {std::tuple<int, size_t, size_t, std::string>{partitioning_scheme::lsh, 8, 5, "equal"},
                                                 {partitioning_scheme::lsh_sim, 257, 8, "skewed"},
                                                 {partitioning_scheme::similarity, 4, 3, "equal"}})
    {
        stress_data const data{approach, n, np, dist};
        try
        {
            data.partition();
            ADD_FAILURE() << "Expected std::runtime_error for approach " << approach;
        }
        catch (std::runtime_error const & error)
        {
            EXPECT_NE(std::string{error.what()}.find("does not fit into any partition"), std::string::npos)
                << error.what();
        }
    }
}

// weighted_fold subtracted the cardinality of the next unassigned large user bin instead of the removed one, counted
// accepted small user bins twice, could remove all user bins of a partition, and its cursors could cross. It crashed
// for valid inputs, e.g. 257 user bins in 4 partitions.
TEST(phibf_regression_test, weighted_fold_crash)
{
    for (auto const & [n, np, dist] : {std::tuple<size_t, size_t, std::string>{257, 4, "skewed"}, {64, 4, "uniform"}})
    {
        stress_data const data{chopper::layout::phibf::partitioning_scheme::weighted_fold, n, np, dist};
        expect_each_user_bin_assigned_once(data.partition(), n);
    }

    // About two user bins per partition: some partitions stay empty, which is reported.
    stress_data const data{chopper::layout::phibf::partitioning_scheme::weighted_fold, 33, 16, "uniform"};
    EXPECT_THROW(data.partition(), std::runtime_error);
}

TEST(phibf_regression_test, weighted_fold)
{
    // 8 user bins, 2 partitions: 125 k-mers and 4 user bins per partition.
    stress_data data{chopper::layout::phibf::partitioning_scheme::weighted_fold, 8, 2, "equal"};
    data.cardinalities = {100, 90, 10, 10, 10, 10, 10, 10};

    // Partition 0 takes the large user bins 100 and 90 (score |1 - 190/125| + |1 - 2/4| = 1.02).
    // Removing 90 and adding small user bins while the score improves: 110/2 (0.62), 120/3 (0.29), 130/4 (0.04);
    // a fourth one would be worse (140/5: 0.37). 0.04 < 1.02, so this is accepted.
    // Removing 100 as well cannot improve on that (at best 30/4: 0.76). The rest goes to partition 1.
    auto const partitions = data.partition();
    ASSERT_EQ(partitions.size(), 2u);

    auto partition_cardinalities = [&data](std::vector<size_t> const & partition)
    {
        std::vector<size_t> result;
        for (size_t const user_bin : partition)
            result.push_back(data.cardinalities[user_bin]);
        std::ranges::sort(result, std::ranges::greater{});
        return result;
    };

    EXPECT_EQ(partition_cardinalities(partitions[0]), (std::vector<size_t>{100, 10, 10, 10}));
    EXPECT_EQ(partition_cardinalities(partitions[1]), (std::vector<size_t>{90, 10, 10, 10}));
    expect_each_user_bin_assigned_once(partitions, 8);
}

namespace
{

// The number of technical bins of the top-level IBF of a layout.
size_t top_level_technical_bins(seqan::hibf::layout::layout const & layout)
{
    size_t result{};
    for (auto const & user_bin : layout.user_bins)
    {
        if (user_bin.previous_TB_indices.empty()) // single or split bin on the top level
            result = std::max(result, user_bin.storage_TB_id + user_bin.number_of_technical_bins);
        else // merged bin on the top level
            result = std::max(result, user_bin.previous_TB_indices[0] + 1);
    }
    return result;
}

std::vector<seqan::hibf::layout::layout> execute_phibf(phibf_data & data, std::filesystem::path const & layout_file)
{
    data.config.output_filename = layout_file;
    std::vector<std::vector<std::string>> filenames(phibf_data::number_of_user_bins, {"ub.fa"});
    chopper::layout::execute(data.config, filenames, data.sketches, data.minHash_sketches);

    std::ifstream layout_stream{layout_file};
    return std::get<2>(chopper::layout::read_layouts_file(layout_stream));
}

} // namespace

// The partitioned HIBF ignored a tmax given by the user and always chose tmax based on each partition's size.
TEST(phibf_execute_test, tmax)
{
    seqan3::test::tmp_directory tmp_dir{};
    size_t const number_of_partitions{3};

    // tmax not set: about 80 user bins per partition, next_multiple_of_64(ceil(sqrt(80))) = 64.
    {
        phibf_data data{number_of_partitions, chopper::layout::phibf::partitioning_scheme::sorted};
        data.config.hibf_config.tmax = 128;
        auto const layouts = execute_phibf(data, tmp_dir.path() / "heuristic.layout");
        ASSERT_EQ(layouts.size(), number_of_partitions);
        for (auto const & layout : layouts)
            EXPECT_EQ(top_level_technical_bins(layout), 64u);
    }

    // tmax set: every partition uses it.
    {
        phibf_data data{number_of_partitions, chopper::layout::phibf::partitioning_scheme::sorted};
        data.config.hibf_config.tmax = 128;
        data.config.tmax_is_set = true;
        auto const layouts = execute_phibf(data, tmp_dir.path() / "tmax.layout");
        ASSERT_EQ(layouts.size(), number_of_partitions);
        for (auto const & layout : layouts)
            EXPECT_EQ(top_level_technical_bins(layout), 128u);
    }
}
