// ---------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
// ---------------------------------------------------------------------------------------------------

#include <algorithm>
#include <cassert>
#include <cinttypes>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include <hibf/contrib/robin_hood.hpp>

#include <chopper/configuration.hpp>
#include <chopper/layout/determine_best_number_of_technical_bins.hpp>
#include <chopper/layout/execute.hpp>
#include <chopper/layout/hibf_statistics.hpp>
#include <chopper/layout/output.hpp>
#include <chopper/next_multiple_of_64.hpp>
#include <chopper/sketch/output.hpp>

#include <hibf/layout/compute_layout.hpp>
#include <hibf/layout/layout.hpp>
#include <hibf/misc/divide_and_ceil.hpp>
#include <hibf/sketch/compute_sketches.hpp>
#include <hibf/sketch/hyperloglog.hpp>
#include <hibf/sketch/toolbox.hpp>

namespace chopper::layout
{

uint64_t lsh_hash_the_sketch(std::vector<uint64_t> const & sketch, size_t const number_of_hashes_to_consider)
{
    // lets just compute the sum
    uint64_t sum{0};

    for (size_t i = 0; i < number_of_hashes_to_consider; ++i)
        sum += sketch[i];

    return sum;
}

// minHash_sketches data structure:
// Vector L1 : number of user bins
// Vector L2 : number_of_max_minHash_sketches (LSH ADD+OR parameter b)
// Vector L3 : minHash_sketche_size (LSH ADD+OR parameter r)
std::vector<std::vector<size_t>> initital_LSH_partitioning(std::vector<std::vector<std::vector<uint64_t>>> const & minHash_sketches,
                         std::vector<size_t> const & cardinalities,
                         size_t const average_technical_bin_size,
                         chopper::configuration const & config)
{
    assert(!minHash_sketches.empty());
    assert(!minHash_sketches[0].empty());
    assert(!minHash_sketches[0][0].empty());

    size_t const number_of_user_bins{cardinalities.size()};
    assert(number_of_user_bins == minHash_sketches.size());
    size_t const number_of_max_minHash_sketches{minHash_sketches[0].size()}; // LSH ADD+OR parameter b
    size_t const minHash_sketche_size{minHash_sketches[0][0].size()};        // LSH ADD+OR parameter r

    size_t number_of_clusters{number_of_user_bins}; // initially, each UB is a separate cluster

    std::vector<size_t> cluster_cardinality(number_of_user_bins);

    // initialise clusters
    // cluster are either
    // 1) of size 1; containing an id != position where the id points to the cluster it has been moved to
    //    e.g. cluster[Y] = {Z} (Y has been moved into Z, so Z could look likes this cluster[Z] = {Z, Y})
    // 2) of size >= 1; with the first entry beging id == position (a valid cluster)
    //    e.g. cluster[X] = {X}       // valid singleton
    //    e.g. cluster[X] = {X, a, b, c, ...}   // valid cluster with more joined entries
    // The clusters could me moved recursively, s.t.
    // cluster[A] = {B}
    // cluster[B] = {C}
    // cluster[C] = {C, A, B} // is valid cluster since cluster[C][0] == C; contains A and B

    std::vector<std::vector<size_t>> clusters(number_of_clusters);
    size_t current_max_cluster_size{0};

    for (size_t user_bin_idx = 0; user_bin_idx < number_of_user_bins; ++user_bin_idx)
    {
        clusters[user_bin_idx] = {user_bin_idx};
        cluster_cardinality[user_bin_idx] = cardinalities[user_bin_idx];
        current_max_cluster_size = std::max(current_max_cluster_size, cardinalities[user_bin_idx]);
    }

    // refine clusters
    size_t current_sketch_index{0};
std::cout << "Start clustering with threshold average_technical_bin_size: " << average_technical_bin_size << std::endl;
    while (current_max_cluster_size < average_technical_bin_size && /*number_of_clusters / static_cast<double>(number_of_user_bins) > 0.5 &&*/
           current_sketch_index < number_of_max_minHash_sketches) // I want to cluster 10%?
    {
std::cout << "Current number of clusters: " << number_of_clusters;

        // fill table
        robin_hood::unordered_flat_map<uint64_t, std::vector<size_t>> table;

        size_t processed_user_bins{0}; // only for sanity check
        for (size_t pos = 0; pos < clusters.size(); ++pos)
        {
            auto const & current = clusters[pos];

            assert(!current.empty());
            size_t const representative_id = current[0];

            // sanity check for cluster
            assert(((current.size() == 1) && (representative_id != pos)) /*current cluster has been moved*/||
                   (representative_id == pos) /*valid cluster, representative id at start equals pos*/);

            if (representative_id != pos) // cluster has been moved somewhere else, don't process
                continue;

            for (size_t const user_bin_idx : current)
            {
                ++processed_user_bins;
                uint64_t const key = lsh_hash_the_sketch(minHash_sketches[user_bin_idx][current_sketch_index], minHash_sketche_size);
                table[key].push_back(representative_id); // insert representative for all user bins
            }
        }
        assert(processed_user_bins == number_of_user_bins); // all user bins should've been processed by one of the clusters

        // read out table
        for (auto & [key, list] : table)
        {
            assert(!list.empty());

            // uniquify list. Since I am inserting representative_idx's into the table, the same number can
            // be inserted into multiple splots, and multiple times in the same slot.
            std::sort(list.begin(), list.end());
            auto const end = std::unique(list.begin(), list.end());
            auto const begin = list.begin();
            size_t current_id = *begin;

            // Now combine all clusters into the first.

            // 1) find the representative cluster
            std::reference_wrapper<std::vector<uint64_t>> representative = clusters[current_id];
            assert(!representative.get().empty());
            assert(((representative.get().size() == 1) && (representative.get()[0] != current_id)) /*representative.get() has been moved*/||
                   (representative.get()[0] == current_id) /*valid cluster,  representative id at start equals pos*/);

            // It can happen, that the representative has already been joined with another cluster
            // e.g.
            // [key1] = {0,11}  // then clusters[11] is merged into clusters[0]
            // [key2] = {11,13} // now I want to merge clusters[13] into clusters[11] but the latter has been moved
            // when cluster[11] has been moved. The entry in cluster[11] is the id of the cluster it has been merged to
            // so we can recursively find a valid cluster
            while (representative.get()[0] != current_id)
            {
                assert(representative.get().size() == 1); // solely contains key where cluster has been moved to
                current_id = representative.get()[0];
                representative = clusters[current_id]; // replace by next cluster
                assert(!representative.get().empty());
            }

            // Now I can be sure that I found the correct representative cluster to merge into
            assert(representative.get()[0] == clusters[representative.get()[0]][0]);
            // Since I requested a user bin id represented by '*begin', at least '*begin' has to be inside the cluster:
            // assert((representative.get().size() == 1) || (std::find(representative.get().begin(), representative.get().end(), *begin) != representative.get().end()));

            for (auto current = begin + 1; current < end; ++current)
            {
                assert((*begin) != (*current));
                // For every other entry in the list, it can happen that I already joined that list with another
                // e.g.
                // [key1] = {0,11}  // then clusters[11] is merged into clusters[0]
                // [key2] = {0, 2, 11} // now I want to do it again
                // I can identify whether it is a moved cluster or not by checking the first entry.
                // e.g. clusters[X] = {X, a, b, c, ...} // valid cluster because clusters[X][0] == X
                // e.g. clusters[Y] = {Z} // moved cluster, because clusters[Y][0] != Y
                current_id = *current;
                std::reference_wrapper<std::vector<uint64_t>> next = clusters[current_id];
                assert(!next.get().empty());
                assert(((next.get().size() == 1) && (next.get()[0] != current_id)) /*next.get() has been moved*/||
                       (next.get()[0] == current_id) /*valid cluster,  next id at start equals pos*/);
                assert(!clusters[next.get()[0]].empty());

                while (next.get()[0] != current_id)
                {
                    assert(next.get().size() == 1); // solely contains key where cluster has been moved to
                    current_id = next.get()[0];
                    next = clusters[current_id];
                    assert(!next.get().empty());
                }

                // Since I requested a user bin id represented by '*current', at least '*current' has to be inside the cluster:
                // assert((next.get().size() == 1) || (std::find(next.get().begin(), next.get().end(), *current) != next.get().end()));

                if (next.get()[0] == representative.get()[0]) // already joined
                    continue;

                // otherwise join
                assert(representative.get()[0] == clusters[representative.get()[0]][0]); // valid cluster
                assert(next.get()[0] == clusters[next.get()[0]][0]); // valid cluster
                representative.get().insert(representative.get().end(), next.get().begin(), next.get().end());
                cluster_cardinality[representative.get()[0]] += cluster_cardinality[next.get().front()];
                assert(representative.get().size() >= 1);
                assert(representative.get()[0] == clusters[representative.get()[0]][0]); // still valid cluster
                next.get() = {representative.get()[0]}; // replace with identifier where cluster has been moved to
                assert(next.get().size() == 1);
                --number_of_clusters;
            }

            current_max_cluster_size = *std::ranges::max_element(cluster_cardinality);
        }

        ++current_sketch_index;

std::cout << " and after this clustering step there are: " << number_of_clusters << "with max cluster size" << current_max_cluster_size << std::endl;
    }

    // clusters are done
    // remove non_valid clusters
    for (size_t pos = 0; pos < clusters.size(); ++pos)
    {
        assert(((clusters[pos].size() == 1) && (clusters[pos][0] != pos)) || (clusters[pos][0] == pos));

        if (clusters[pos][0] != pos)
            clusters[pos].clear();

        // sort ub wihtin a cluster by cardinality
        std::sort(clusters[pos].begin(), clusters[pos].end(), [&cardinalities](auto const & v1, auto const & v2){ return cardinalities[v2] < cardinalities[v1]; });
    }

    // push largest p clusters to the front
    auto cluster_size_cmp = [&cardinalities](auto const & v1, auto const & v2)
    {
        if (v2.size() == v1.size() && !v2.empty())
            return  cardinalities[v2[0]] < cardinalities[v1[0]];
        return v2.size() < v1.size();
    };
    std::partial_sort(clusters.begin(), clusters.begin() + config.number_of_partitions, clusters.end(), cluster_size_cmp);

    // after filling up the partitions with the biggest clusters, sort the clusters by cardinality of the biggest ub
    // s.t. that euqally sizes ub are assigned after each other and the small stuff is added at last.
    // the largest ub is already at the start because of former sorting.
    auto compare_cardinality_and_move_empty_clusters_to_the_end = [&cardinalities](auto const & v1, auto const & v2)
    {
        if (v1.empty())
            return false; // v1 can never be larger than v2 then

        if (v2.empty()) // and v1 is not, since the first if would catch
            return true;

        return cardinalities[v2[0]] < cardinalities[v1[0]];
    };
    std::sort(clusters.begin() + config.number_of_partitions, clusters.end(), compare_cardinality_and_move_empty_clusters_to_the_end);

    assert(clusters[0].size() >= clusters[1].size()); // sanity check
    assert(cardinalities[clusters[config.number_of_partitions][0]] >= cardinalities[clusters[config.number_of_partitions + 1][0]]); // sanity check

// debug
    for (size_t cidx = 1; cidx < clusters.size(); ++cidx)
    {
        if (clusters[cidx - 1].empty() && !clusters[cidx].empty()) // once empty - always empty; all empty clusters should be at the end
            throw std::runtime_error{"sorting did not work"};
    }
// debug

    return clusters;
}


bool find_best_partition(chopper::configuration const & config,
                         size_t const cardinality_per_part,
                         std::vector<size_t> const & cluster,
                         std::vector<size_t> const & cardinalities,
                         std::vector<seqan::hibf::sketch::hyperloglog> const & sketches,
                         std::vector<std::vector<size_t>> & positions,
                         std::vector<seqan::hibf::sketch::hyperloglog> & partition_sketches,
                         std::vector<size_t> & partition_cardinality)
{
    seqan::hibf::sketch::hyperloglog const current_sketch = [&sketches, &cluster, &config]()
    {
        seqan::hibf::sketch::hyperloglog result{config.hibf_config.sketch_bits};

        for (size_t const user_bin_idx : cluster)
            result.merge(sketches[user_bin_idx]);

        return result;
    }();

    size_t const current_cardinality = [&cardinalities, &cluster]()
    {
        size_t result{0};

        for (size_t const user_bin_idx : cluster)
            result += cardinalities[user_bin_idx];

        return result;
    }();

    // TODO: afterwads check if I should again merge by how little the effective text ratio grows

    // search best partition fit by similarity
    // similarity here is defined as:
    // "whose (<-partition) effective text size is subsumed most by the current user bin"
    // or in other words:
    // "which partition has the largest intersection with user bin b compared to its own (partition) size."
    double best_subsume_ratio{0.0};
    size_t best_p{0};
    bool best_p_found{false};

    for (size_t p = 0; p < config.number_of_partitions; ++p)
    {
        seqan::hibf::sketch::hyperloglog tmp = current_sketch;
        tmp.merge(partition_sketches[p]);
        size_t const tmp_estimate = tmp.estimate();
        size_t const current_partition_size = partition_sketches[p].estimate();
        assert(tmp_estimate >= current_partition_size);
        size_t const change = tmp_estimate - current_partition_size;
        size_t const intersection = current_sketch.estimate() - change;
        double const subsume_ratio = static_cast<double>(intersection) / current_partition_size;

        if (subsume_ratio > best_subsume_ratio &&
            partition_cardinality[p] < cardinality_per_part &&
            (current_cardinality + partition_cardinality[p]) < (cardinality_per_part * 1.2) )
        {
            best_subsume_ratio = subsume_ratio;
            best_p = p;
            best_p_found = true;
        }
    }

    if (!best_p_found)
        return false;

    // now that we know which partition fits best (`best_p`), add those indices to it
    for (size_t const user_bin_idx : cluster)
    {
        positions[best_p].push_back(user_bin_idx);
        partition_cardinality[best_p] += cardinalities[user_bin_idx];
    }
    partition_sketches[best_p].merge(current_sketch);

    return true;
}

void partition_user_bins(chopper::configuration const & config,
                         std::vector<size_t> const & cardinalities,
                         std::vector<seqan::hibf::sketch::hyperloglog> const & sketches,
                         std::vector<std::vector<std::vector<uint64_t>>> const & minHash_sketches,
                         std::vector<std::vector<size_t>> & partitions)
{
    // all approaches need sorted positions
    std::vector<size_t> const sorted_positions = [&cardinalities]()
    {
        std::vector<size_t> ps;
        ps.resize(cardinalities.size());
        std::iota(ps.begin(), ps.end(), 0);
        seqan::hibf::sketch::toolbox::sort_by_cardinalities(cardinalities, ps);
        return ps;
    }();

    size_t const sum_of_cardinalities = [&cardinalities]()
    {
        size_t sum{0};
        for (size_t const card : cardinalities)
            sum += card;
        return sum;
    }();

    if (config.partitioning_approach == partitioning_scheme::blocked)
    {
        size_t const u_bins_per_part = seqan::hibf::divide_and_ceil(cardinalities.size(), config.number_of_partitions);
        size_t const block_size =
            std::min(u_bins_per_part,
                     chopper::next_multiple_of_64(static_cast<uint16_t>(std::ceil(std::sqrt(u_bins_per_part)))));

        size_t current_part{0u};
        size_t current_block_count{0};

        for (size_t const current_user_bin_id : sorted_positions)
        {
            partitions[current_part].push_back(current_user_bin_id);
            ++current_block_count;

            if (current_block_count >= block_size)
            {
                current_block_count = 0;
                ++current_part;
                if (current_part == config.number_of_partitions) // we need to circle back to the first partition
                    current_part = 0;
            }
        }
    }
    else if (config.partitioning_approach == partitioning_scheme::sorted)
    {
        size_t const cardinality_per_part =
            seqan::hibf::divide_and_ceil(sum_of_cardinalities, config.number_of_partitions);

        size_t current_cardinality{0u};
        size_t current_part{0};

        for (size_t const current_user_bin_id : sorted_positions)
        {
            partitions[current_part].push_back(current_user_bin_id);
            current_cardinality += cardinalities[current_user_bin_id];

            if (current_cardinality >= cardinality_per_part)
            {
                current_cardinality = 0;
                ++current_part;
            }
        }
    }
    else if (config.partitioning_approach == partitioning_scheme::folded)
    {
        size_t const cardinality_per_part_halved =
            seqan::hibf::divide_and_ceil(sum_of_cardinalities, config.number_of_partitions * 2);

        size_t current_cardinality{0u};
        std::vector<size_t> const parts = [&config]()
        {
            size_t const len{config.number_of_partitions};
            std::vector<size_t> result(len * 2);
            std::iota(result.begin(), result.begin() + len, 0);
            std::copy(result.rbegin() + len, result.rend(), result.begin() + len);
            return result;
        }();
        size_t current_part{0};

        for (size_t const current_user_bin_id : sorted_positions)
        {
            partitions[parts[current_part]].push_back(current_user_bin_id);
            current_cardinality += cardinalities[current_user_bin_id];

            if (current_cardinality >= cardinality_per_part_halved)
            {
                current_cardinality = 0;
                ++current_part;
            }
        }
    }
    else if (config.partitioning_approach == partitioning_scheme::weighted_fold)
    {
        size_t const cardinality_per_part =
            seqan::hibf::divide_and_ceil(sum_of_cardinalities, config.number_of_partitions);
        size_t const u_bins_per_part = seqan::hibf::divide_and_ceil(cardinalities.size(), config.number_of_partitions);

        size_t current_big_pos{0};                          // the next largest user bin to assign to a partition
        size_t current_small_pos{cardinalities.size() - 1}; // the next small user bin

        for (size_t current_part = 0; current_part + 1 < config.number_of_partitions; ++current_part)
        {
            size_t current_cardinality{0};
            std::vector<size_t> small_bins;
            size_t new_small_bin_addition{0};

            auto compute_score = [&]()
            {
                double const weight = static_cast<double>(current_cardinality) / cardinality_per_part;
                double const amount =
                    static_cast<double>(partitions[current_part].size() + small_bins.size() + new_small_bin_addition)
                    / u_bins_per_part;
                return std::abs(1.0 - weight) + std::abs(1.0 - amount);
            };

            // first add all large bins that fit
            while (current_cardinality < cardinality_per_part)
            {
                partitions[current_part].push_back(sorted_positions[current_big_pos]);
                current_cardinality += cardinalities[sorted_positions[current_big_pos]];
                ++current_big_pos;
            }

            double local_optimum = compute_score();

            // then remove big bins and add small bins until a local optima is reached
            while (true)
            {
                size_t const cache_last_small_pos{current_small_pos};
                // remove a big user bin and fill the partition with small user bins
                current_cardinality -= cardinalities[sorted_positions[current_big_pos]];

                // can we further improve the ratio by adding more small bins?
                double improved_score{};
                do
                {
                    improved_score = compute_score();
                    current_cardinality += cardinalities[sorted_positions[current_small_pos]];
                    --current_small_pos;
                    ++new_small_bin_addition;
                }
                while (compute_score() < improved_score); // smaller is better
                // remove overstep
                ++current_small_pos;
                current_cardinality -= cardinalities[sorted_positions[current_small_pos]];
                --new_small_bin_addition;

                if (local_optimum < compute_score()) // score would increase. Stop
                {
                    current_small_pos = cache_last_small_pos;
                    break;
                }
                else // update
                {
                    partitions[current_part].pop_back();
                    --current_big_pos;
                    for (size_t pos = cache_last_small_pos; pos > current_small_pos; --pos)
                        small_bins.push_back(sorted_positions[pos]);
                }
            }
            partitions[current_part].insert(partitions[current_part].end(), small_bins.begin(), small_bins.end());
        }

        // remaining user bins go to last partition
        while (current_big_pos <= current_small_pos)
        {
            partitions[config.number_of_partitions - 1].push_back(sorted_positions[current_big_pos]);
            ++current_big_pos;
        }
    }
    else if (config.partitioning_approach == partitioning_scheme::similarity)
    {
        uint8_t const sketch_bits{config.hibf_config.sketch_bits};
        std::vector<seqan::hibf::sketch::hyperloglog> partition_sketches(config.number_of_partitions,
                                                                         seqan::hibf::sketch::hyperloglog(sketch_bits));
        std::vector<size_t> partition_cardinality(config.number_of_partitions, 0u);

        size_t const cardinality_per_part =
            seqan::hibf::divide_and_ceil(sum_of_cardinalities, config.number_of_partitions);
        size_t const u_bins_per_part = seqan::hibf::divide_and_ceil(cardinalities.size(), config.number_of_partitions);
        size_t const block_size =
            std::min(u_bins_per_part,
                     chopper::next_multiple_of_64(static_cast<uint16_t>(std::ceil(std::sqrt(u_bins_per_part)))));
        size_t const number_of_blocks = seqan::hibf::divide_and_ceil(cardinalities.size(), block_size);

        // don't move from largest to smallest but pick the next block to process randomly.
        // this probably leads to more evenly distributed partitions (evenly in terms of number of user bins)
        std::vector<size_t> indices(number_of_blocks);
        std::iota(indices.begin(), indices.end(), 0);
        std::random_device shuffle_random_device;
        std::mt19937 shuffle_engine(shuffle_random_device());
        std::shuffle(indices.begin(), indices.end(), shuffle_engine);

        // initialise partitions with the first random config.number_of_partitions blocks
        assert(number_of_blocks >= config.number_of_partitions);
        for (size_t p = 0; p < config.number_of_partitions; ++p)
        {
            size_t const i = indices[p];
            size_t const end = (i == (indices.size() - 1) && (cardinalities.size() % block_size != 0)) ? (cardinalities.size() % block_size) : block_size;
            for (size_t x = 0; x < end; ++x)
            {
                assert(block_size * i + x < sorted_positions.size());
                partition_sketches[p].merge(sketches[sorted_positions[block_size * i + x]]);
                partition_cardinality[p] += cardinalities[sorted_positions[block_size * i + x]];
                partitions[p].push_back(sorted_positions[block_size * i + x]);
            }
        }

        auto has_been_processed_in_init = [&config,&indices, &block_size](size_t const i)
        {
            bool result{false};
            for (size_t p = 0; p < config.number_of_partitions; ++p)
                result |= ((indices[p] * block_size) == i);
            return result;
        };

        // assign the rest by similarity, user in by user bin (single bins, no blocks)
        for (size_t i = 0; i < sorted_positions.size(); ++i)
        {
            if (has_been_processed_in_init(i))
            {
                i += (block_size - 1); // -1 because there will be an increment after continue
                continue;
            }

            find_best_partition(config,
                                cardinality_per_part,
                                {sorted_positions[i]},
                                cardinalities,
                                sketches,
                                partitions,
                                partition_sketches,
                                partition_cardinality);
        }
    }
    else if (config.partitioning_approach == partitioning_scheme::lsh)
    {
        uint8_t const sketch_bits{config.hibf_config.sketch_bits};
        std::vector<seqan::hibf::sketch::hyperloglog> partition_sketches(config.number_of_partitions,
                                                                         seqan::hibf::sketch::hyperloglog(sketch_bits));
        std::vector<size_t> partition_cardinality(config.number_of_partitions, 0u);

        size_t const cardinality_per_part =
            seqan::hibf::divide_and_ceil(sum_of_cardinalities, config.number_of_partitions);

        // initial partitioning using locality sensitive hashing (LSH)
        std::vector<std::vector<size_t>> const clusters = initital_LSH_partitioning(minHash_sketches, cardinalities, cardinality_per_part, config);

        // initialise partitions with the first p largest clusters (initital_LSH_partitioning sorts by size)
        size_t cidx{0}; // current cluster index
        for (size_t p = 0; p < config.number_of_partitions; ++p)
        {
            assert(!clusters[cidx].empty());
            bool p_has_been_incremented{false};

            for (size_t uidx = 0; uidx < clusters[cidx].size(); ++uidx)
            {
                size_t const user_bin_idx = clusters[cidx][uidx];
                // if a single cluster already exceeds the cardinality_per_part,
                // then the remaining user bins of the cluster must spill over into the next partition
                if (partition_cardinality[p] > cardinality_per_part)
                {
                    ++p;
                    p_has_been_incremented = true;
                    assert(p < config.number_of_partitions);
                }

                partition_sketches[p].merge(sketches[user_bin_idx]);
                partition_cardinality[p] += cardinalities[user_bin_idx];
                partitions[p].push_back(user_bin_idx);
            }

            if (p_has_been_incremented && partition_cardinality[p] < cardinality_per_part)
            {
                // p will be incremeted again in the next iteration of this for loop
                // since we have incremented p previously because of some user bins that spilled over in this currnet p
                // we want to decrease p here, s.t. we stay at the same p.
                // this ensures that the current p doesn't only have a small rest that just didn't fit into the
                //  previous partition
                assert(p > 0);
                --p;
            }

            ++cidx;
        }

        // assign the rest by similarity
        assert(config.hibf_config.number_of_user_bins == clusters.size());
        for (; cidx < config.hibf_config.number_of_user_bins; ++cidx)
        {
            auto const & cluster = clusters[cidx];

            if (cluster.empty()) // non valid clusters have been removed. Since list was sorted, we can abort
                break;

            bool const cluster_could_be_added = find_best_partition(config, cardinality_per_part, cluster, cardinalities, sketches, partitions, partition_sketches, partition_cardinality);

            // if the cluster is too large to be added. Add single user bins
            if (!cluster_could_be_added)
            {
                for (size_t const user_bin_idx : cluster)
                {
                    bool const found = find_best_partition(config, cardinality_per_part, {user_bin_idx}, cardinalities, sketches, partitions, partition_sketches, partition_cardinality);
                    assert(found); // there should always be at least one partition that has enough space left
                }
            }
        }
    }
    else if (config.partitioning_approach == partitioning_scheme::lsh_sim)
    {
        uint8_t const sketch_bits{config.hibf_config.sketch_bits};
        std::vector<seqan::hibf::sketch::hyperloglog> partition_sketches(config.number_of_partitions,
                                                                         seqan::hibf::sketch::hyperloglog(sketch_bits));
        std::vector<size_t> partition_cardinality(config.number_of_partitions, 0u);

        size_t const cardinality_per_part =
            seqan::hibf::divide_and_ceil(sum_of_cardinalities, config.number_of_partitions);

        // initial partitioning using locality sensitive hashing (LSH)
        std::vector<std::vector<size_t>> const clusters = initital_LSH_partitioning(minHash_sketches, cardinalities, cardinality_per_part, config);

        // initialise partitions with the first p largest clusters (initital_LSH_partitioning sorts by size)
        size_t cidx{0}; // current cluster index
        for (size_t p = 0; p < config.number_of_partitions; ++p)
        {
            assert(!clusters[cidx].empty());
            bool p_has_been_incremented{false};

            for (size_t uidx = 0; uidx < clusters[cidx].size(); ++uidx)
            {
                size_t const user_bin_idx = clusters[cidx][uidx];
                // if a single cluster already exceeds the cardinality_per_part,
                // then the remaining user bins of the cluster must spill over into the next partition
                if (partition_cardinality[p] > cardinality_per_part)
                {
                    ++p;
                    p_has_been_incremented = true;
                    assert(p < config.number_of_partitions);
                }

                partition_sketches[p].merge(sketches[user_bin_idx]);
                partition_cardinality[p] += cardinalities[user_bin_idx];
                partitions[p].push_back(user_bin_idx);
            }

            if (p_has_been_incremented && partition_cardinality[p] < cardinality_per_part)
            {
                // p will be incremeted again in the next iteration of this for loop
                // since we have incremented p previously because of some user bins that spilled over in this currnet p
                // we want to decrease p here, s.t. we stay at the same p.
                // this ensures that the current p doesn't only have a small rest that just didn't fit into the
                //  previous partition
                assert(p > 0);
                --p;
            }

            ++cidx;
        }

        // assign the rest, ub by ub, by similarity
        assert(config.hibf_config.number_of_user_bins == clusters.size());
        std::vector<size_t> remaining_ubs{};
        for (; cidx < config.hibf_config.number_of_user_bins; ++cidx)
        {
            auto const & cluster = clusters[cidx];

            if (cluster.empty()) // non valid clusters have been removed. Since list was sorted, we can abort
                break;

            remaining_ubs.insert(remaining_ubs.end(), cluster.begin(), cluster.end());
        }

        std::sort(remaining_ubs.begin(), remaining_ubs.end(),
            [&cardinalities] (size_t const ub1, size_t const ub2) { return cardinalities[ub2] < cardinalities[ub1]; });

        assert(remaining_ubs.size() < 2 || cardinalities[remaining_ubs[0]] >= cardinalities[remaining_ubs[1]]);

        for (size_t const user_bin_idx : remaining_ubs)
        {
            bool const found = find_best_partition(config, cardinality_per_part, {user_bin_idx}, cardinalities, sketches, partitions, partition_sketches, partition_cardinality);
            assert(found); // there should always be at least one partition that has enough space left
        }
    }

    // sanity check:
    size_t sum{0};
    for (auto const & p : partitions)
        sum += p.size();

    if (sum != sketches.size())
    {
        std::string str{"Not all user bins have been assigned to the "};
        str += std::to_string(partitions.size());
        str += " partitions! (";
        str += std::to_string(sum);
        str += "/";
        str += std::to_string(sketches.size());
        str += ")\n";
        for (auto const & p : partitions)
        {
            str += "[";
            for (auto const h : p)
            {
                str += std::to_string(h);
                str += ",";
            }
            str.back() = ']';
            str += '\n';
        }

        throw std::logic_error{str};
    }
}

int execute(chopper::configuration & config, std::vector<std::vector<std::string>> const & filenames)
{
    assert(config.hibf_config.number_of_user_bins > 0);

    if (config.hibf_config.disable_estimate_union)
        config.hibf_config.disable_rearrangement = true;

    if (config.hibf_config.tmax == 0) // no tmax was set by the user on the command line
    {
        // Set default as sqrt(#samples). Experiments showed that this is a reasonable default.
        if (size_t number_samples = config.hibf_config.number_of_user_bins;
            number_samples >= 1ULL << 32) // sqrt is bigger than uint16_t
            throw std::invalid_argument{"Too many samples. Please set a tmax (see help via `-hh`)."}; // GCOVR_EXCL_LINE
        else
            config.hibf_config.tmax =
                chopper::next_multiple_of_64(static_cast<uint16_t>(std::ceil(std::sqrt(number_samples))));
    }
    else if (config.hibf_config.tmax % 64 != 0)
    {
        config.hibf_config.tmax = chopper::next_multiple_of_64(config.hibf_config.tmax);
        std::cerr << "[CHOPPER LAYOUT WARNING]: Your requested number of technical bins was not a multiple of 64. "
                  << "Due to the architecture of the HIBF, it will use up space equal to the next multiple of 64 "
                  << "anyway, so we increased your number of technical bins to " << config.hibf_config.tmax << ".\n";
    }

    if (config.number_of_partitions < 2) // 0 == unset == single HIBF, 1 == single HIBF
    {
        seqan::hibf::layout::layout hibf_layout;
        std::vector<seqan::hibf::sketch::hyperloglog> sketches;

        seqan::hibf::concurrent_timer compute_sketches_timer{};
        seqan::hibf::concurrent_timer union_estimation_timer{};
        seqan::hibf::concurrent_timer rearrangement_timer{};
        seqan::hibf::concurrent_timer dp_algorithm_timer{};

        if (config.determine_best_tmax)
        {
            std::tie(hibf_layout, sketches) = determine_best_number_of_technical_bins(config);
        }
        else
        {
            std::vector<size_t> kmer_counts;

            compute_sketches_timer.start();
            seqan::hibf::sketch::compute_sketches(config.hibf_config, kmer_counts, sketches);
            compute_sketches_timer.stop();

            std::vector<size_t> positions = [&kmer_counts]()
            {
                std::vector<size_t> ps;
                ps.resize(kmer_counts.size());
                std::iota(ps.begin(), ps.end(), 0);
                return ps;
            }(); // GCOVR_EXCL_LINE

            dp_algorithm_timer.start();
            hibf_layout = seqan::hibf::layout::compute_layout(config.hibf_config,
                                                              kmer_counts,
                                                              sketches,
                                                              std::move(positions),
                                                              union_estimation_timer,
                                                              rearrangement_timer);
            dp_algorithm_timer.stop();

            if (config.output_verbose_statistics)
            {
                size_t dummy{};
                chopper::layout::hibf_statistics global_stats{config, sketches, kmer_counts};
                global_stats.hibf_layout = hibf_layout;
                global_stats.print_header_to(std::cout);
                global_stats.print_summary_to(dummy, std::cout);
            }
        }

        if (!config.disable_sketch_output)
        {
            if (!std::filesystem::exists(config.sketch_directory))
                std::filesystem::create_directory(config.sketch_directory);

            assert(filenames.size() == sketches.size());
            for (size_t i = 0; i < filenames.size(); ++i)
                sketch::write_sketch_file(filenames[i][0], sketches[i], config);
        }

        // brief Write the output to the layout file.
        std::ofstream fout{config.output_filename};
        chopper::layout::write_user_bins_to(filenames, fout);
        config.write_to(fout);
        hibf_layout.write_to(fout);

        if (!config.output_timings.empty())
        {
            std::ofstream output_stream{config.output_timings};
            output_stream << std::fixed << std::setprecision(2);
            output_stream << "sketching_in_seconds\t"
                          << "layouting_in_seconds\t"
                          << "union_estimation_in_seconds\t"
                          << "rearrangement_in_seconds\n";
            output_stream << compute_sketches_timer.in_seconds() << '\t';
            output_stream << dp_algorithm_timer.in_seconds() << '\t';
            output_stream << union_estimation_timer.in_seconds() << '\t';
            output_stream << rearrangement_timer.in_seconds() << '\t';
        }
    }
    else
    {
        std::vector<size_t> cardinalities;
        std::vector<seqan::hibf::sketch::hyperloglog> sketches;
        std::vector<std::vector<size_t>> positions(config.number_of_partitions); // asign positions for each partition
        std::vector<std::vector<std::vector<uint64_t>>> minHash_sketches;

        // compute sketches of all user bins
        seqan::hibf::concurrent_timer compute_sketches_timer{};
        compute_sketches_timer.start();
        seqan::hibf::sketch::compute_sketches_with_minhash(config.hibf_config, cardinalities, sketches, minHash_sketches);
        compute_sketches_timer.stop();

        partition_user_bins(config, cardinalities, sketches, minHash_sketches, positions);

        std::vector<seqan::hibf::layout::layout> hibf_layouts(config.number_of_partitions); // multiple layouts

#pragma omp parallel for schedule(dynamic) num_threads(config.hibf_config.threads)
        for (size_t i = 0; i < config.number_of_partitions; ++i)
        {
            seqan::hibf::concurrent_timer union_estimation_timer{};
            seqan::hibf::concurrent_timer rearrangement_timer{};
            seqan::hibf::concurrent_timer dp_algorithm_timer{};

            // reset tmax to fit number of user bins in layout
            auto local_hibf_config = config.hibf_config; // every thread needs to set individual tmax
            local_hibf_config.tmax =
                chopper::next_multiple_of_64(static_cast<uint16_t>(std::ceil(std::sqrt(positions[i].size()))));

            dp_algorithm_timer.start();
            hibf_layouts[i] = seqan::hibf::layout::compute_layout(local_hibf_config,
                                                                  cardinalities,
                                                                  sketches,
                                                                  std::move(positions[i]),
                                                                  union_estimation_timer,
                                                                  rearrangement_timer);
            dp_algorithm_timer.stop();

            if (!config.output_timings.empty())
            {
                std::ofstream output_stream{config.output_timings, std::ios_base::app};
                output_stream << std::fixed << std::setprecision(2);
                output_stream << "sketching_in_seconds\t"
                              << "layouting_in_seconds\t"
                              << "union_estimation_in_seconds\t"
                              << "rearrangement_in_seconds\n";
                output_stream << compute_sketches_timer.in_seconds() << '\t';
                output_stream << dp_algorithm_timer.in_seconds() << '\t';
                output_stream << union_estimation_timer.in_seconds() << '\t';
                output_stream << rearrangement_timer.in_seconds() << '\t';
            }
        }

        // brief Write the output to the layout file.
        std::ofstream fout{config.output_filename};
        chopper::layout::write_user_bins_to(filenames, fout);
        config.write_to(fout);

        for (size_t i = 0; i < config.number_of_partitions; ++i)
            hibf_layouts[i].write_to(fout);
    }

    return 0;
}

} // namespace chopper::layout
