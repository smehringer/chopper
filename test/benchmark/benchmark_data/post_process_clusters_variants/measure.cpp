// ---------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
// ---------------------------------------------------------------------------------------------------

// Compares two variants of post_process_clusters (src/layout/partition_user_bins.cpp) on one synthetic input.
// Variant A sorts the first tmax clusters by size (as in chopper), variant B the first number_of_remaining_tbs.
// Must be linked with a copy of partition_user_bins.cpp that has variant_b.patch applied. See README.md.
//
// Usage: measure <n> <family_divisor> <sigma> <seed> <median> <min> <max> <threads>
// Prints one tab-separated line; run.sh prints the header.

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <map>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include <chopper/configuration.hpp>
#include <chopper/layout/fast_layout.hpp>
#include <chopper/layout/hibf_statistics.hpp>
#include <chopper/layout/partition_user_bins.hpp>

#include <hibf/layout/compute_fpr_correction.hpp>
#include <hibf/layout/compute_relaxed_fpr_correction.hpp>
#include <hibf/layout/layout.hpp>
#include <hibf/sketch/compute_sketches.hpp>
#include <hibf/sketch/estimate_kmer_counts.hpp>

// Defined in the patched partition_user_bins.cpp.
extern bool chopper_variant_b;

namespace
{

// splitmix64 finaliser: turns consecutive integers into well-distributed hashes.
uint64_t scramble(uint64_t x)
{
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

// A user bin consists of a random subset of its family's k-mer pool plus unique k-mers.
struct user_bin_spec
{
    uint64_t family;
    uint64_t pool_size;
    uint32_t keep_threshold; // Keep a pool k-mer if a hash of (user bin, k-mer) & 0xffff is below this.
    uint64_t unique;
};

struct top_metrics
{
    size_t split_tbs{};     // Technical bins holding a user bin that spans at least two technical bins.
    double max_corrected{}; // Largest FPR-corrected technical bin.
    std::vector<std::vector<size_t>> partitions{};
};

// Partitions the top level and computes its metrics.
top_metrics top_level(chopper::configuration const & config,
                      std::vector<size_t> const & positions,
                      std::vector<size_t> const & cardinalities,
                      std::vector<seqan::hibf::sketch::hyperloglog> const & sketches,
                      std::vector<seqan::hibf::sketch::minhashes> const & minhashes)
{
    top_metrics m{};
    m.partitions.resize(config.hibf_config.tmax);
    chopper::layout::partition_user_bins(config, positions, cardinalities, sketches, minhashes, m.partitions);

    auto const split_corr =
        seqan::hibf::layout::compute_fpr_correction({.fpr = config.hibf_config.maximum_fpr,
                                                     .hash_count = config.hibf_config.number_of_hash_functions,
                                                     .t_max = config.hibf_config.tmax});
    double const relaxed = seqan::hibf::layout::compute_relaxed_fpr_correction(
        {.fpr = config.hibf_config.maximum_fpr,
         .relaxed_fpr = config.hibf_config.relaxed_fpr,
         .hash_count = config.hibf_config.number_of_hash_functions});

    std::map<size_t, size_t> occurrences{};
    for (auto const & p : m.partitions)
        if (p.size() == 1)
            ++occurrences[p[0]];

    for (auto const & p : m.partitions)
    {
        if (p.empty())
            continue;
        if (p.size() > 1)
        {
            seqan::hibf::sketch::hyperloglog u{config.hibf_config.sketch_bits};
            for (size_t ub : p)
                u.merge(sketches[ub]);
            m.max_corrected = std::max(m.max_corrected, u.estimate() * relaxed);
        }
        else
        {
            size_t const k = occurrences[p[0]];
            if (k > 1)
                ++m.split_tbs;
            m.max_corrected = std::max(m.max_corrected, cardinalities[p[0]] * split_corr[k] / k);
        }
    }
    return m;
}

} // namespace

int main(int argc, char ** argv)
{
    if (argc != 9)
    {
        std::fprintf(stderr, "Usage: %s <n> <family_divisor> <sigma> <seed> <median> <min> <max> <threads>\n", argv[0]);
        return 1;
    }
    size_t const n = std::stoul(argv[1]);
    size_t const family_divisor = std::stoul(argv[2]);
    double const sigma = std::stod(argv[3]);
    uint64_t const seed = std::stoul(argv[4]);
    double const median = std::stod(argv[5]);
    uint64_t const min_pool = std::stoul(argv[6]);
    uint64_t const max_pool = std::stoul(argv[7]);
    size_t const threads = std::stoul(argv[8]);

    // Input generation. Family pool sizes are log-normal; each user bin keeps 60 % to 98 % of its family's pool.
    size_t const families = std::max<size_t>(1, n / family_divisor);
    std::mt19937_64 rng{seed * 1000003 + n * 31 + family_divisor * 7 + static_cast<uint64_t>(sigma * 10)};
    std::lognormal_distribution<double> size_dist{std::log(median), sigma};
    std::uniform_real_distribution<double> sim_dist{0.6, 0.98};

    std::vector<uint64_t> family_pool(families);
    for (auto & b : family_pool)
        b = std::clamp<uint64_t>(static_cast<uint64_t>(size_dist(rng)), min_pool, max_pool);

    std::vector<user_bin_spec> specs(n);
    for (size_t u = 0; u < n; ++u)
    {
        uint64_t const f = rng() % families;
        double const sim = sim_dist(rng);
        specs[u] = {f,
                    family_pool[f],
                    static_cast<uint32_t>(sim * 65536),
                    static_cast<uint64_t>((1.0 - sim) * family_pool[f]) + 500};
    }

    chopper::configuration config{};
    config.hibf_config.number_of_user_bins = n;
    config.hibf_config.threads = threads;
    config.hibf_config.input_fn = [&specs, seed](size_t const u, seqan::hibf::insert_iterator it)
    {
        auto const & sp = specs[u];
        for (uint64_t j = 0; j < sp.pool_size; ++j)
            if ((scramble(seed ^ (static_cast<uint64_t>(u) << 32) ^ j) & 0xffff) < sp.keep_threshold)
                it = scramble((sp.family << 36) + j + (seed << 60));
        for (uint64_t i = 0; i < sp.unique; ++i)
            it = scramble((1ULL << 63) | (static_cast<uint64_t>(u) << 34) | i) ^ seed;
    };
    config.hibf_config.validate_and_set_defaults(); // tmax as chosen by the command line interface

    auto const t0 = std::chrono::steady_clock::now();
    std::vector<seqan::hibf::sketch::hyperloglog> sketches{};
    std::vector<seqan::hibf::sketch::minhashes> minhashes{};
    seqan::hibf::sketch::compute_sketches(config.hibf_config, sketches, minhashes);
    std::vector<size_t> cardinalities{};
    seqan::hibf::sketch::estimate_kmer_counts(sketches, cardinalities);
    std::vector<size_t> positions(n);
    std::iota(positions.begin(), positions.end(), 0u);
    double const t_sketch = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();

    // Both variants on the same sketches.
    top_metrics top[2];
    size_t size[2];
    double cost[2];
    double t_layout[2];
    size_t depth[2];
    for (int v = 0; v < 2; ++v)
    {
        chopper_variant_b = (v == 1);
        top[v] = top_level(config, positions, cardinalities, sketches, minhashes);

        auto const t1 = std::chrono::steady_clock::now();
        seqan::hibf::layout::layout hibf_layout{};
        chopper::layout::fast_layout(config, positions, cardinalities, sketches, minhashes, hibf_layout);
        t_layout[v] = std::chrono::duration<double>(std::chrono::steady_clock::now() - t1).count();

        depth[v] = 0;
        for (auto const & ub : hibf_layout.user_bins)
            depth[v] = std::max(depth[v], ub.previous_TB_indices.size() + 1);

        chopper::layout::hibf_statistics stats{config, sketches, cardinalities};
        stats.hibf_layout = hibf_layout;
        size[v] = stats.total_hibf_size_in_byte();
        cost[v] = stats.expected_HIBF_query_cost;
    }

    std::printf("%zu\t%zu\t%.1f\t%lu\t%.0f\t%lu\t%lu\t" // n families sigma seed median min max
                "%zu\t%zu\t%zu\t%d\t%.0f\t%.0f\t"       // tmax s_A s_B same_top topmax_A topmax_B
                "%zu\t%zu\t%.4f\t%.4f\t%zu\t%zu\t"      // size_A size_B cost_A cost_B depth_A depth_B
                "%.1f\t%.1f\t%.1f\n",                   // t_sketch t_A t_B
                n,
                families,
                sigma,
                static_cast<unsigned long>(seed),
                median,
                static_cast<unsigned long>(min_pool),
                static_cast<unsigned long>(max_pool),
                config.hibf_config.tmax,
                top[0].split_tbs,
                top[1].split_tbs,
                top[0].partitions == top[1].partitions ? 1 : 0,
                top[0].max_corrected,
                top[1].max_corrected,
                size[0],
                size[1],
                cost[0],
                cost[1],
                depth[0],
                depth[1],
                t_sketch,
                t_layout[0],
                t_layout[1]);
}
