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
