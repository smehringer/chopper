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
#include <optional>

#include <hibf/contrib/robin_hood.hpp>

#include <chopper/configuration.hpp>
#include <chopper/layout/determine_best_number_of_technical_bins.hpp>
#include <chopper/layout/execute.hpp>
#include <chopper/layout/hibf_statistics.hpp>
#include <chopper/layout/output.hpp>
#include <chopper/next_multiple_of_64.hpp>
#include <chopper/sketch/output.hpp>
#include <chopper/lsh.hpp>

#include <hibf/layout/compute_layout.hpp>
#include <hibf/layout/layout.hpp>
#include <hibf/misc/divide_and_ceil.hpp>
#include <hibf/misc/iota_vector.hpp>
#include <hibf/sketch/estimate_kmer_counts.hpp> // for estimate_kmer_counts
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

auto LSH_fill_hashtable(std::vector<Cluster> const & clusters,
                        std::vector<seqan::hibf::sketch::minhashes> const & minHash_sketches,
                        size_t const current_sketch_index,
                        size_t const current_number_of_sketch_hashes)
{
    robin_hood::unordered_flat_map<uint64_t, std::vector<size_t>> table;

    size_t processed_user_bins{0}; // only for sanity check
    for (size_t pos = 0; pos < clusters.size(); ++pos)
    {
        auto const & current = clusters[pos];
        assert(current.is_valid(pos));

        if (current.has_been_moved()) // cluster has been moved somewhere else, don't process
            continue;

        for (size_t const user_bin_idx : current.contained_user_bins())
        {
            ++processed_user_bins;
            uint64_t const key = lsh_hash_the_sketch(minHash_sketches[user_bin_idx].table[current_sketch_index], current_number_of_sketch_hashes);
            table[key].push_back(current.id()); // insert representative for all user bins
        }
    }
    assert(processed_user_bins == clusters.size()); // all user bins should've been processed by one of the clusters

    return table;
}

auto LSH_fill_hashtable(std::vector<MultiCluster> const & clusters,
                        std::vector<seqan::hibf::sketch::minhashes> const & minHash_sketches,
                        size_t const current_sketch_index,
                        size_t const current_number_of_sketch_hashes)
{
    robin_hood::unordered_flat_map<uint64_t, std::vector<size_t>> table;

    size_t processed_user_bins{0}; // only for sanity check
    for (size_t pos = 0; pos < clusters.size(); ++pos)
    {
        auto const & current = clusters[pos];
        assert(current.is_valid(pos));

        if (current.has_been_moved()) // cluster has been moved somewhere else, don't process
            continue;

        for (auto const & similarity_cluster : current.contained_user_bins())
        {
            for (size_t const user_bin_idx : similarity_cluster)
            {
                ++processed_user_bins;
                uint64_t const key = lsh_hash_the_sketch(minHash_sketches[user_bin_idx].table[current_sketch_index], current_number_of_sketch_hashes);
                table[key].push_back(current.id()); // insert representative for all user bins
            }
        }
    }
    assert(processed_user_bins == clusters.size()); // all user bins should've been processed by one of the clusters

    return table;
}

// minHash_sketches data structure:
// Vector L1 : number of user bins
// Vector L2 : number_of_max_minHash_sketches (LSH ADD+OR parameter b)
// Vector L3 : minHash_sketche_size (LSH ADD+OR parameter r)
std::vector<Cluster> initital_LSH_partitioning(std::vector<seqan::hibf::sketch::minhashes> const & minHash_sketches,
                         std::vector<size_t> const & positions,
                         std::vector<size_t> const & cardinalities,
                         std::vector<seqan::hibf::sketch::hyperloglog> const & sketches,
                         size_t const average_technical_bin_size,
                         chopper::configuration const & config)
{
    assert(!minHash_sketches.empty());
    assert(!minHash_sketches[0].table.empty());
    assert(!minHash_sketches[0].table[0].empty());

    size_t const number_of_user_bins{positions.size()};
    size_t const number_of_max_minHash_sketches{seqan::hibf::sketch::minhashes::num_sketches}; // LSH ADD+OR parameter b
    size_t const minHash_sketche_size{seqan::hibf::sketch::minhashes::sketch_size};        // LSH ADD+OR parameter r
    seqan::hibf::sketch::hyperloglog const empty_sketch{config.hibf_config.sketch_bits};

    // initialise clusters with a signle user bin per cluster.
    // clusters are either
    // 1) of size 1; containing an id != position where the id points to the cluster it has been moved to
    //    e.g. cluster[Y] = {Z} (Y has been moved into Z, so Z could look likes this cluster[Z] = {Z, Y})
    // 2) of size >= 1; with the first entry beging id == position (a valid cluster)
    //    e.g. cluster[X] = {X}       // valid singleton
    //    e.g. cluster[X] = {X, a, b, c, ...}   // valid cluster with more joined entries
    // The clusters could me moved recursively, s.t.
    // cluster[A] = {B}
    // cluster[B] = {C}
    // cluster[C] = {C, A, B} // is valid cluster since cluster[C][0] == C; contains A and B
    std::vector<Cluster> clusters;
    clusters.reserve(number_of_user_bins);

    std::vector<size_t> current_cluster_cardinality(number_of_user_bins);
    std::vector<seqan::hibf::sketch::hyperloglog> current_cluster_sketches(number_of_user_bins, empty_sketch);
    size_t current_max_cluster_size{0};
    size_t current_number_of_sketch_hashes{minHash_sketche_size}; // start with high r but decrease it iteratively
    size_t current_sketch_index{0};
    size_t current_number_of_clusters{number_of_user_bins}; // initially, each UB is a separate cluster

    for (size_t pos = 0; pos < number_of_user_bins; ++pos)
    {
        clusters.emplace_back(pos, positions[pos]); // id, user_bin_idx
        current_cluster_cardinality[pos] = cardinalities[positions[pos]];
        current_cluster_sketches[pos] = sketches[positions[pos]];
        current_max_cluster_size = std::max(current_max_cluster_size, cardinalities[positions[pos]]);
    }

    // refine clusters
//std::cout << "Start clustering with threshold average_technical_bin_size: " << average_technical_bin_size << std::endl;
    while (current_max_cluster_size < average_technical_bin_size && /*number_of_clusters / static_cast<double>(number_of_user_bins) > 0.5 &&*/
           current_sketch_index < number_of_max_minHash_sketches) // I want to cluster 10%?
    {
//std::cout << "Current number of clusters: " << current_number_of_clusters;

        // fill LSH collision hashtable
        robin_hood::unordered_flat_map<uint64_t, std::vector<size_t>> table =
            LSH_fill_hashtable(clusters, minHash_sketches, current_sketch_index, current_number_of_sketch_hashes);

        // read out LSH collision hashtable
        // for each present key, if the list contains more than one cluster, we merge everything contained in the list
        // into the first cluster, since those clusters collide and should be joined in the same bucket
        for (auto & [key, list] : table)
        {
            assert(!list.empty());

            // uniquify list. Since I am inserting representative_idx's into the table, the same number can
            // be inserted into multiple splots, and multiple times in the same slot.
            std::sort(list.begin(), list.end());
            auto const end = std::unique(list.begin(), list.end());
            auto const begin = list.begin();

            if (end - begin <= 1) // nothing to do here
                continue;

            // Now combine all clusters into the first.

            // 1) find the representative cluster to merge everything else into
            // It can happen, that the representative has already been joined with another cluster
            // e.g.
            // [key1] = {0,11}  // then clusters[11] is merged into clusters[0]
            // [key2] = {11,13} // now I want to merge clusters[13] into clusters[11] but the latter has been moved
            size_t const representative_cluster_id = LSH_find_representative_cluster(clusters, *begin);
            auto & representative_cluster = clusters[representative_cluster_id];
            assert(representative_cluster.is_valid(representative_cluster_id));
            assert(representative_cluster.id() == clusters[representative_cluster.id()].id());

            for (auto current = begin + 1; current < end; ++current)
            {
                // For every other entry in the list, it can happen that I already joined that list with another
                // e.g.
                // [key1] = {0,11}  // then clusters[11] is merged into clusters[0]
                // [key2] = {0, 2, 11} // now I want to do it again
                size_t const next_cluster_id = LSH_find_representative_cluster(clusters, *current);
                auto & next_cluster = clusters[next_cluster_id];

                if (next_cluster.id() == representative_cluster.id()) // already joined
                    continue;

                next_cluster.move_to(representative_cluster); // otherwise join next_cluster into representative_cluster
                assert(next_cluster.size() == 0);
                assert(next_cluster.has_been_moved());
                assert(representative_cluster.size() > 1); // there should be at least two user bins now
                assert(representative_cluster.is_valid(representative_cluster_id)); // and it should still be valid

                current_cluster_sketches[representative_cluster.id()].merge(current_cluster_sketches[next_cluster.id()]);
                current_cluster_cardinality[representative_cluster.id()] = current_cluster_sketches[representative_cluster.id()].estimate();

                --current_number_of_clusters;
            }

            current_max_cluster_size = *std::ranges::max_element(current_cluster_cardinality);
        }

        ++current_sketch_index;

//std::cout << " and after this clustering step there are: " << current_number_of_clusters << "with max cluster size" << current_max_cluster_size << std::endl;
    }

    return clusters;
}

std::vector<Cluster> very_similar_LSH_partitioning(std::vector<seqan::hibf::sketch::minhashes> const & minHash_sketches,
                         std::vector<size_t> const & positions,
                         std::vector<size_t> const & cardinalities,
                         std::vector<seqan::hibf::sketch::hyperloglog> const & sketches,
                         size_t const average_technical_bin_size,
                         chopper::configuration const & config)
{
    assert(!minHash_sketches.empty());
    assert(!minHash_sketches[0].table.empty());
    assert(!minHash_sketches[0].table[0].empty());

    size_t const number_of_user_bins{positions.size()};
    assert(number_of_user_bins == minHash_sketches.size());
    size_t const number_of_max_minHash_sketches{3};                     // LSH ADD+OR parameter b
    size_t const minHash_sketche_size{minHash_sketches[0].table[0].size()};   // LSH ADD+OR parameter r
    seqan::hibf::sketch::hyperloglog const empty_sketch{config.hibf_config.sketch_bits};

    // initialise clusters with a signle user bin per cluster.
    // clusters are either
    // 1) of size 1; containing an id != position where the id points to the cluster it has been moved to
    //    e.g. cluster[Y] = {Z} (Y has been moved into Z, so Z could look likes this cluster[Z] = {Z, Y})
    // 2) of size >= 1; with the first entry beging id == position (a valid cluster)
    //    e.g. cluster[X] = {X}       // valid singleton
    //    e.g. cluster[X] = {X, a, b, c, ...}   // valid cluster with more joined entries
    // The clusters could me moved recursively, s.t.
    // cluster[A] = {B}
    // cluster[B] = {C}
    // cluster[C] = {C, A, B} // is valid cluster since cluster[C][0] == C; contains A and B
    std::vector<Cluster> clusters;
    clusters.reserve(number_of_user_bins);

    std::vector<size_t> current_cluster_cardinality(number_of_user_bins);
    std::vector<seqan::hibf::sketch::hyperloglog> current_cluster_sketches(number_of_user_bins, empty_sketch);
    size_t current_max_cluster_size{0};
    size_t current_number_of_sketch_hashes{minHash_sketche_size}; // start with high r but decrease it iteratively
    size_t current_sketch_index{0};
    size_t current_number_of_clusters{number_of_user_bins}; // initially, each UB is a separate cluster

    for (size_t pos = 0; pos < number_of_user_bins; ++pos)
    {
        clusters.emplace_back(pos, positions[pos]);
        current_cluster_cardinality[pos] = cardinalities[positions[pos]];
        current_cluster_sketches[pos] = sketches[positions[pos]];
        current_max_cluster_size = std::max(current_max_cluster_size, cardinalities[positions[pos]]);
    }

    // refine clusters
//std::cout << "Start clustering with threshold average_technical_bin_size: " << average_technical_bin_size << std::endl;
    while (current_max_cluster_size < average_technical_bin_size && /*number_of_clusters / static_cast<double>(number_of_user_bins) > 0.5 &&*/
           current_sketch_index < number_of_max_minHash_sketches) // I want to cluster 10%?
    {
//std::cout << "Current number of clusters: " << current_number_of_clusters;

        // fill LSH collision hashtable
        robin_hood::unordered_flat_map<uint64_t, std::vector<size_t>> table =
            LSH_fill_hashtable(clusters, minHash_sketches, current_sketch_index, current_number_of_sketch_hashes);

        // read out LSH collision hashtable
        // for each present key, if the list contains more than one cluster, we merge everything contained in the list
        // into the first cluster, since those clusters collide and should be joined in the same bucket
        for (auto & [key, list] : table)
        {
            assert(!list.empty());

            // uniquify list. Since I am inserting representative_idx's into the table, the same number can
            // be inserted into multiple splots, and multiple times in the same slot.
            std::sort(list.begin(), list.end());
            auto const end = std::unique(list.begin(), list.end());
            auto const begin = list.begin();

            if (end - begin <= 1) // nothing to do here
                continue;

            // Now combine all clusters into the first.

            // 1) find the representative cluster to merge everything else into
            // It can happen, that the representative has already been joined with another cluster
            // e.g.
            // [key1] = {0,11}  // then clusters[11] is merged into clusters[0]
            // [key2] = {11,13} // now I want to merge clusters[13] into clusters[11] but the latter has been moved
            size_t const representative_cluster_id = LSH_find_representative_cluster(clusters, *begin);
            auto & representative_cluster = clusters[representative_cluster_id];
            assert(representative_cluster.is_valid(representative_cluster_id));
            assert(representative_cluster.id() == clusters[representative_cluster.id()].id());

            for (auto current = begin + 1; current < end; ++current)
            {
                // For every other entry in the list, it can happen that I already joined that list with another
                // e.g.
                // [key1] = {0,11}  // then clusters[11] is merged into clusters[0]
                // [key2] = {0, 2, 11} // now I want to do it again
                size_t const next_cluster_id = LSH_find_representative_cluster(clusters, *current);
                auto & next_cluster = clusters[next_cluster_id];

                if (next_cluster.id() == representative_cluster.id()) // already joined
                    continue;

                next_cluster.move_to(representative_cluster); // otherwise join next_cluster into representative_cluster
                assert(next_cluster.size() == 0);
                assert(next_cluster.has_been_moved());
                assert(representative_cluster.size() > 1); // there should be at least two user bins now
                assert(representative_cluster.is_valid(representative_cluster_id)); // and it should still be valid

                current_cluster_sketches[representative_cluster.id()].merge(current_cluster_sketches[next_cluster.id()]);
                current_cluster_cardinality[representative_cluster.id()] = current_cluster_sketches[representative_cluster.id()].estimate();

                --current_number_of_clusters;
            }

            current_max_cluster_size = *std::ranges::max_element(current_cluster_cardinality);
        }

        ++current_sketch_index;

//std::cout << " and after this clustering step there are: " << current_number_of_clusters << "with max cluster size" << current_max_cluster_size << std::endl;
    }

    return clusters;
}

std::vector<MultiCluster> most_distant_LSH_partitioning(std::vector<Cluster> const & initial_clusters,
                         std::vector<seqan::hibf::sketch::minhashes> const & minHash_sketches,
                         chopper::configuration const & config)
{
    assert(!minHash_sketches.empty());
    assert(!minHash_sketches[0].table.empty());
    assert(!minHash_sketches[0].table[0].empty());

    size_t const number_of_user_bins{initial_clusters.size()};
    size_t const number_of_max_minHash_sketches{minHash_sketches[0].table.size()}; // LSH ADD+OR parameter b
    // size_t const minHash_sketche_size{minHash_sketches[0][0].size()};   // LSH ADD+OR parameter r

    size_t current_number_of_sketch_hashes{5};
    size_t current_sketch_index{0};

    size_t current_number_of_clusters = [&initial_clusters] ()
    {
        size_t number{0};
        for (auto const & cluster : initial_clusters)
        {
            if (cluster.size())
                ++number;
        }
        return number;
    }();

    std::vector<MultiCluster> clusters;
    clusters.reserve(number_of_user_bins);

    for (size_t pos = 0; pos < number_of_user_bins; ++pos)
    {
        clusters.emplace_back(initial_clusters[pos]);
        assert(clusters[pos].is_valid(pos));
    }

    // refine clusters. I want to cluster as much as possible to initialise partitions with very distant clusters
    while (current_number_of_clusters > (config.hibf_config.tmax * 2) &&
           current_sketch_index < number_of_max_minHash_sketches)
    {
//std::cout << "[Dist] Current number of clusters: " << current_number_of_clusters << " sketch size:" << current_number_of_sketch_hashes << " ";

        // fill LSH collision hashtable
        robin_hood::unordered_flat_map<uint64_t, std::vector<size_t>> table =
            LSH_fill_hashtable(clusters, minHash_sketches, current_sketch_index, current_number_of_sketch_hashes);

        // read out LSH collision hashtable
        // for each present key, if the list contains more than one cluster, we merge everything contained in the list
        // into the first cluster, since those clusters collide and should be joined in the same bucket
        for (auto & [key, list] : table)
        {
            assert(!list.empty());

            // uniquify list. Since I am inserting representative_idx's into the table, the same number can
            // be inserted into multiple splots, and multiple times in the same slot.
            std::sort(list.begin(), list.end());
            auto const end = std::unique(list.begin(), list.end());
            auto const begin = list.begin();

            if (end - begin <= 1) // nothing to do here
                continue;

            // Now combine all clusters into the first.

            // 1) find the representative cluster to merge everything else into
            // It can happen, that the representative has already been joined with another cluster
            // e.g.
            // [key1] = {0,11}  // then clusters[11] is merged into clusters[0]
            // [key2] = {11,13} // now I want to merge clusters[13] into clusters[11] but the latter has been moved
            size_t const representative_cluster_id = LSH_find_representative_cluster(clusters, *begin);
            auto & representative_cluster = clusters[representative_cluster_id];
            assert(representative_cluster.is_valid(representative_cluster_id));
            assert(representative_cluster.id() == clusters[representative_cluster.id()].id());

            for (auto current = begin + 1; current < end; ++current)
            {
                // For every other entry in the list, it can happen that I already joined that list with another
                // e.g.
                // [key1] = {0,11}  // then clusters[11] is merged into clusters[0]
                // [key2] = {0, 2, 11} // now I want to do it again
                size_t const next_cluster_id = LSH_find_representative_cluster(clusters, *current);
                auto & next_cluster = clusters[next_cluster_id];

                if (next_cluster.id() == representative_cluster.id()) // already joined
                    continue;

                next_cluster.move_to(representative_cluster); // otherwise join next_cluster into representative_cluster
                assert(next_cluster.size() == 0);
                assert(next_cluster.has_been_moved());
                assert(representative_cluster.size() > 1); // there should be at least two user bins now
                assert(representative_cluster.is_valid(representative_cluster_id)); // and it should still be valid

                --current_number_of_clusters;
            }
        }

        ++current_sketch_index;

        if (current_sketch_index % 4 == 0 && current_number_of_sketch_hashes > 1)
            --current_number_of_sketch_hashes;

//std::cout << " and after this clustering step there are: " << current_number_of_clusters << std::endl;
    }

    return clusters;
}

void post_process_clusters(std::vector<Cluster> & clusters,
                           std::vector<size_t> const & cardinalities,
                           chopper::configuration const & config)
{
    // clusters are done. Start post processing
    // since post processing involves re-ordering the clusters, the moved_to_cluster_id value of a cluster will not
    // refer to the position of the cluster in the `clusters` vecto anymore but the cluster with the resprive id()
    // would neet to be found

    for (size_t pos = 0; pos < clusters.size(); ++pos)
    {
        assert(clusters[pos].is_valid(pos));
        clusters[pos].sort_by_cardinality(cardinalities);
    }

    // push largest p clusters to the front
    auto cluster_size_cmp = [&cardinalities](auto const & v1, auto const & v2)
    {
        if (v2.size() == v1.size() && !v2.empty())
            return  cardinalities[v2.contained_user_bins()[0]] < cardinalities[v1.contained_user_bins()[0]];
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

        return cardinalities[v2.contained_user_bins()[0]] < cardinalities[v1.contained_user_bins()[0]];
    };
    std::sort(clusters.begin() + config.number_of_partitions, clusters.end(), compare_cardinality_and_move_empty_clusters_to_the_end);

    assert(clusters[0].size() >= clusters[1].size()); // sanity check
    assert(cardinalities[clusters[config.number_of_partitions].contained_user_bins()[0]] >= cardinalities[clusters[config.number_of_partitions + 1].contained_user_bins()[0]]); // sanity check

// debug
    for (size_t cidx = 1; cidx < clusters.size(); ++cidx)
    {
        if (clusters[cidx - 1].empty() && !clusters[cidx].empty()) // once empty - always empty; all empty clusters should be at the end
            throw std::runtime_error{"sorting did not work"};
    }
// debug
}

void post_process_clusters(std::vector<MultiCluster> & clusters,
                           std::vector<size_t> const & cardinalities,
                           chopper::configuration const & config)
{
    // clusters are done. Start post processing
    // since post processing involves re-ordering the clusters, the moved_to_cluster_id value of a cluster will not
    // refer to the position of the cluster in the `clusters` vecto anymore but the cluster with the resprive id()
    // would neet to be found

    for (size_t pos = 0; pos < clusters.size(); ++pos)
    {
        assert(clusters[pos].is_valid(pos));
        clusters[pos].sort_by_cardinality(cardinalities);
    }

    // push largest p clusters to the front
    auto cluster_size_cmp = [&cardinalities](auto const & mc1, auto const & mc2)
    {
        if (mc1.empty())
            return false; // mc1 can never be larger than mc2 then
        if (mc2.empty()) // and mc1 is not, since the first if would catch
            return true;

        auto const & largest_c1 = mc1.contained_user_bins()[0]; // first cluster is the largest because of sorting above
        auto const & largest_c2 = mc2.contained_user_bins()[0]; // first cluster is the largest because of sorting above

        assert(!largest_c1.empty()); // should never be empty
        assert(!largest_c2.empty()); // should never be empty

        if (largest_c2.size() == largest_c1.size())
            return  cardinalities[largest_c2[0]] < cardinalities[largest_c1[0]];
        return largest_c2.size() < largest_c1.size();
    };
    std::partial_sort(clusters.begin(), clusters.begin() + config.number_of_partitions, clusters.end(), cluster_size_cmp);

    auto move_empty_clusters_to_the_end = [](auto const & v1, auto const & v2)
    {
        return v2.size() < v1.size();
    };
    std::sort(clusters.begin() + config.number_of_partitions, clusters.end(), move_empty_clusters_to_the_end);

// debug sanity check
    for (size_t cidx = 1; cidx < clusters.size(); ++cidx)
    {
        if (clusters[cidx - 1].empty() && !clusters[cidx].empty()) // once empty - always empty; all empty clusters should be at the end
            throw std::runtime_error{"sorting did not work"};
    }
// debug
}

std::vector<Cluster> LSH_partitioning(std::vector<seqan::hibf::sketch::minhashes> const & minHash_sketches,
                                      std::vector<size_t> const & positions,
                                      std::vector<size_t> const & cardinalities,
                                      std::vector<seqan::hibf::sketch::hyperloglog> const & sketches,
                                      size_t const average_technical_bin_size,
                                      chopper::configuration const & config)
{
    std::vector<Cluster> clusters = initital_LSH_partitioning(minHash_sketches, positions, cardinalities, sketches, average_technical_bin_size, config);
    post_process_clusters(clusters, cardinalities, config);
    return clusters;
}

std::vector<MultiCluster> sim_dist_LSH_partitioning(std::vector<seqan::hibf::sketch::minhashes> const & minHash_sketches,
                                      std::vector<size_t> const & positions,
                                      std::vector<size_t> const & cardinalities,
                                      std::vector<seqan::hibf::sketch::hyperloglog> const & sketches,
                                      size_t const average_technical_bin_size,
                                      chopper::configuration const & config)
{
    std::vector<Cluster> clusters = very_similar_LSH_partitioning(minHash_sketches, positions, cardinalities, sketches, average_technical_bin_size, config);
    std::vector<MultiCluster> multi_clusters = most_distant_LSH_partitioning(clusters, minHash_sketches, config);
    post_process_clusters(multi_clusters, cardinalities, config);
    return multi_clusters;
}

bool find_best_partition(chopper::configuration const & config,
                         size_t & corrected_estimate_per_part,
                         std::vector<size_t> const & cluster,
                         std::vector<size_t> const & cardinalities,
                         std::vector<seqan::hibf::sketch::hyperloglog> const & sketches,
                         std::vector<std::vector<size_t>> & positions,
                         std::vector<seqan::hibf::sketch::hyperloglog> & partition_sketches,
                         std::vector<size_t> & max_partition_cardinality,
                         std::vector<size_t> & min_partition_cardinality)
{
    seqan::hibf::sketch::hyperloglog const current_sketch = [&sketches, &cluster, &config]()
    {
        seqan::hibf::sketch::hyperloglog result{config.hibf_config.sketch_bits};

        for (size_t const user_bin_idx : cluster)
            result.merge(sketches[user_bin_idx]);

        return result;
    }();

    size_t const max_card = [&cardinalities,  &cluster] ()
    {
        size_t max{0};

        for (size_t const user_bin_idx : cluster)
            max = std::max(max, cardinalities[user_bin_idx]);

        return max;
    }();

    // TODO: afterwads check if I should again merge by how little the effective text ratio grows

    // search best partition fit by similarity
    // similarity here is defined as:
    // "whose (<-partition) effective text size is subsumed most by the current user bin"
    // or in other words:
    // "which partition has the largest intersection with user bin b compared to its own (partition) size."
    // double best_subsume_ratio{0.0};
    size_t smallest_change{std::numeric_limits<size_t>::max()};
    size_t best_p{0};
    bool best_p_found{false};

    auto penalty_lower_level = [&] (size_t const additional_number_of_user_bins, size_t const p)
    {
        size_t min = std::min(max_card, min_partition_cardinality[p]);
        size_t max = std::max(max_card, max_partition_cardinality[p]);

        if (positions[p].size() > config.hibf_config.tmax) // already a third level
        {
            // if there must already be another lower level because the current merged bin contains more than tmax
            // user bins, then the current user bin is very likely stored multiple times. Therefore, the penalty is set
            // to the cardinality of the current user bin times the number of levels, e.g. the number of times this user
            // bin needs to be stored additionally
            size_t const num_ubs_in_merged_bin{positions[p].size() + additional_number_of_user_bins};
            double const levels = std::log(num_ubs_in_merged_bin) / std::log(config.hibf_config.tmax);
            return  static_cast<size_t>(max_card * levels);
        }
        else if (positions[p].size() + additional_number_of_user_bins > config.hibf_config.tmax) // now a third level
        {
            // if the current merged bin contains exactly tmax UBS, adding otherone must
            // result in another lower level. Most likely, the smallest user bin will end up on the lower level
            // therefore the penalty is set to 'min * tmax'
            // of course, there could also be a third level with a lower number of user bins, but this is hard to
            // estimate.
            size_t const penalty = min * config.hibf_config.tmax;
            return penalty;
        }
        else // positions[p].size() + additional_number_of_user_bins < tmax
        {
            // if the new user bin is smaller than all other already contained user bins
            // the waste of space is high if stored in a single technical bin
            if (max_card < min)
                return (max - max_card);
            // if the new user bin is bigger than all other already contained user bins, the IBF size increases
            else if (max_card > max)
                return (max_card - max) * config.hibf_config.tmax;
            // else, end if-else-block and zero is returned
        }

        return (size_t)0u;
    };

    for (size_t p = 0; p < config.number_of_partitions; ++p)
    {
        seqan::hibf::sketch::hyperloglog union_sketch = current_sketch;
        union_sketch.merge(partition_sketches[p]);
        size_t const union_estimate = union_sketch.estimate();
        size_t const current_partition_size = partition_sketches[p].estimate();

        assert(union_estimate >= current_partition_size);
        size_t const penalty_current_bin = union_estimate - current_partition_size;
        size_t const penalty_current_ibf = config.hibf_config.tmax * ((union_estimate <= corrected_estimate_per_part) ? 0u : union_estimate - corrected_estimate_per_part);
        size_t const change = penalty_current_bin + penalty_current_ibf + penalty_lower_level(cluster.size(), p);

        // size_t const intersection = current_sketch.estimate() - change;
        // double const subsume_ratio = static_cast<double>(intersection) / current_partition_size;
//std::cout << "p:" << p << " p-#UBs" << positions[p].size() << " penalty:" <<  penalty(cluster.size(), p) << " change:" << change << " union-current_p:" << (union_estimate - current_partition_size) << " union:" << union_estimate << " current_p:" << current_partition_size << " t:" << corrected_estimate_per_part << std::endl;
        if (change == 0 || /* If there is no penalty at all, this is a best fit even if the partition is "full"*/
            (smallest_change > change /*&& subsume_ratio > best_subsume_ratio &&*/
            /*current_partition_size < corrected_estimate_per_part*/))
        {
//std::cout << "smaller!" << std::endl;
            // best_subsume_ratio = subsume_ratio;
            smallest_change = change;
            best_p = p;
            best_p_found = true;
        }
    }

    if (!best_p_found)
        throw "currently there are no safety measures if a partition is not found because it is very unlikely";

//std::cout << "best_p:" << best_p << std::endl<< std::endl;

    // now that we know which partition fits best (`best_p`), add those indices to it
    for (size_t const user_bin_idx : cluster)
    {
        positions[best_p].push_back(user_bin_idx);
        max_partition_cardinality[best_p] = std::max(max_partition_cardinality[best_p], cardinalities[user_bin_idx]);
        min_partition_cardinality[best_p] = std::min(min_partition_cardinality[best_p], cardinalities[user_bin_idx]);
    }
    partition_sketches[best_p].merge(current_sketch);
    corrected_estimate_per_part = std::max(corrected_estimate_per_part, static_cast<size_t>(partition_sketches[best_p].estimate()));

    return true;
}

bool find_best_sorted_partition(chopper::configuration const & config,
                         size_t & corrected_estimate_per_part,
                         std::vector<size_t> const & cluster,
                         std::vector<size_t> const & cardinalities,
                         std::vector<seqan::hibf::sketch::hyperloglog> const & /*sketches*/,
                         std::vector<std::vector<size_t>> & positions,
                         std::vector<size_t> & partition_sizes,
                         std::vector<size_t> & max_partition_cardinality,
                         std::vector<size_t> & min_partition_cardinality)
{
    size_t max_card = [&cardinalities,  &cluster] ()
    {
        size_t max{0};

        for (size_t const user_bin_idx : cluster)
            max = std::max(max, cardinalities[user_bin_idx]);

        return max;
    }();

    size_t smallest_change{std::numeric_limits<size_t>::max()};
    size_t best_p{0};
    bool best_p_found{false};

    auto penalty_lower_level = [&] (size_t const additional_number_of_user_bins, size_t const p)
    {
        size_t min = std::min(max_card, min_partition_cardinality[p]);
        size_t max = std::max(max_card, max_partition_cardinality[p]);

        if (positions[p].size() > config.hibf_config.tmax) // already a third level
        {
            // if there must already be another lower level because the current merged bin contains more than tmax
            // user bins, then the current user bin is very likely stored multiple times. Therefore, the penalty is set
            // to the cardinality of the current user bin times the number of levels, e.g. the number of times this user
            // bin needs to be stored additionally
            size_t const num_ubs_in_merged_bin{positions[p].size() + additional_number_of_user_bins};
            double const levels = std::log(num_ubs_in_merged_bin) / std::log(config.hibf_config.tmax);
            return  static_cast<size_t>(max_card * levels);
        }
        else if (positions[p].size() + additional_number_of_user_bins > config.hibf_config.tmax)
        {
            // if the current merged bin contains exactly tmax UBS, adding otherone must
            // result in another lower level. Most likely, the smallest user bin will end up on the lower level
            // therefore the penalty is set to 'min * tmax'
            // of course, there could also be a third level with a lower number of user bins, but this is hard to
            // estimate.
            size_t const penalty = min * config.hibf_config.tmax;
            return penalty;
        }
        else // positions[p].size() + additional_number_of_user_bins <= tmax
        {
            // if the new user bin is smaller than all other already contained user bins
            // the waste of space is high if stored in a single technical bin
            if (max_card < min)
                return (max - max_card);
            // if the new user bin is bigger than all other already contained user bins, the IBF size increases
            else if (max_card > max)
                return (max_card - max) * config.hibf_config.tmax;
            // else, end if-else-block and zero is returned
        }

        return (size_t)0u;
    };

    for (size_t p = 0; p < config.number_of_partitions; ++p)
    {
        size_t const change = penalty_lower_level(cluster.size(), p);

        if (change == 0 || /* If there is no penalty at all, this is a best fit even if the partition is "full"*/
            (smallest_change > change && /*subsume_ratio > best_subsume_ratio &&*/
            partition_sizes[p] < corrected_estimate_per_part))
        {
            smallest_change = change;
            best_p = p;
            best_p_found = true;
        }
    }

    if (!best_p_found)
        return false;

//std::cout << "best_p:" << best_p << std::endl<< std::endl;

    // now that we know which partition fits best (`best_p`), add those indices to it
    for (size_t const user_bin_idx : cluster)
    {
        positions[best_p].push_back(user_bin_idx);
        partition_sizes[best_p] += cardinalities[user_bin_idx];
        max_partition_cardinality[best_p] = std::max(max_partition_cardinality[best_p], cardinalities[user_bin_idx]);
        min_partition_cardinality[best_p] = std::min(min_partition_cardinality[best_p], cardinalities[user_bin_idx]);
    }
    corrected_estimate_per_part = std::max(corrected_estimate_per_part, static_cast<size_t>(partition_sizes[best_p]));

    return true;
}

void partition_user_bins(chopper::configuration const & config,
                         std::vector<size_t> const & positions,
                         std::vector<size_t> const & cardinalities,
                         std::vector<seqan::hibf::sketch::hyperloglog> const & sketches,
                         std::vector<seqan::hibf::sketch::minhashes> const & minHash_sketches,
                         std::vector<std::vector<size_t>> & partitions)
{
    // all approaches need sorted positions
    std::vector<size_t> const sorted_positions = [&positions, &cardinalities]()
    {
        std::vector<size_t> ps(positions.begin(), positions.end());
        seqan::hibf::sketch::toolbox::sort_by_cardinalities(cardinalities, ps);
        return ps;
    }();

    size_t const sum_of_cardinalities = [&positions, &cardinalities]()
    {
        size_t sum{0};
        for (size_t const pos : positions)
            sum += cardinalities[pos];
        return sum;
    }();

    size_t const joint_estimate = [&config, &positions, &sketches]()
    {
        seqan::hibf::sketch::hyperloglog sketch{config.hibf_config.sketch_bits};
        for (size_t const pos : positions)
            sketch.merge(sketches[pos]);
        return sketch.estimate();
    }();

    size_t const max_cardinality = [&positions, &cardinalities]()
    {
        size_t max{0};
        for (size_t const pos : positions)
            max = std::max(max, cardinalities[pos]);
        return max;
    }();

    // If the effective text size is very low, it can happen that the joint_estimate divided by the number of partitions
    // is lower than the largest single user bin. But of course, we can never reach a smaller max technical bin size
    // then that of the largest user user bin. Thus we can correct the estimate_per_part beforehand.
    // This way we make sure there is at least 1 LSH clustering step.
    size_t const estimate_per_part = std::max(seqan::hibf::divide_and_ceil(joint_estimate, config.number_of_partitions),
                                              max_cardinality + 1);

//std::cout << "sum_of_cardinalities:" << sum_of_cardinalities << " joint_estimate:" << joint_estimate << std::endl;

    if (config.partitioning_approach == partitioning_scheme::blocked)
    {
        size_t const u_bins_per_part = seqan::hibf::divide_and_ceil(sorted_positions.size(), config.number_of_partitions);
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
        size_t current_tb_size_threshold = estimate_per_part * (1.0 + 2 * static_cast<double>(joint_estimate)/sum_of_cardinalities);

        size_t current_part{0};
        seqan::hibf::sketch::hyperloglog current_sketch{config.hibf_config.sketch_bits};
        std::vector<seqan::hibf::sketch::hyperloglog> partition_sketches(config.number_of_partitions, current_sketch);

        size_t idx{0};
        for (;idx < sorted_positions.size(); ++idx)
        {
            size_t const current_user_bin_id{sorted_positions[idx]};
            partitions[current_part].push_back(current_user_bin_id);
            current_sketch.merge(sketches[current_user_bin_id]);

            if (current_sketch.estimate() >= current_tb_size_threshold)
            {
                partition_sketches[current_part] = current_sketch;
                current_sketch = seqan::hibf::sketch::hyperloglog{config.hibf_config.sketch_bits};
                ++current_part;

                if (current_part >= config.number_of_partitions)
                {
                    current_tb_size_threshold *= 1 + static_cast<double>(idx) / sorted_positions.size();
                    current_part = 0; // fill up from the start agagin as the threshold increased
                }
            }
        }
    }
    else if (config.partitioning_approach == partitioning_scheme::folded)
    {
        size_t const estimate_per_part_halved = seqan::hibf::divide_and_ceil(estimate_per_part, 2u);

        size_t current_tb_size_threshold = estimate_per_part_halved * (1.0 + static_cast<double>(joint_estimate)/sum_of_cardinalities);

        std::vector<size_t> const parts = [&config]()
        {
            size_t const len{config.number_of_partitions};
            std::vector<size_t> result(len * 2);
            std::iota(result.begin(), result.begin() + len, 0);
            std::copy(result.rbegin() + len, result.rend(), result.begin() + len);
            return result;
        }();

        seqan::hibf::sketch::hyperloglog current_sketch{config.hibf_config.sketch_bits};
        size_t current_part{0};
        std::vector<seqan::hibf::sketch::hyperloglog> partition_sketches(config.number_of_partitions, current_sketch);

        size_t idx{0};
        for (;idx < sorted_positions.size(); ++idx)
        {
            partitions[parts[current_part]].push_back(sorted_positions[idx]);
            current_sketch.merge(sketches[sorted_positions[idx]]);

            if (current_sketch.estimate() >= current_tb_size_threshold)
            {
                partition_sketches[parts[current_part]].merge(current_sketch);
                current_sketch = seqan::hibf::sketch::hyperloglog{config.hibf_config.sketch_bits};
                ++current_part;

                if (current_part >= parts.size())
                {
                    current_tb_size_threshold *= 1 + static_cast<double>(idx) / sorted_positions.size();
                    current_part = 0; // fill up from the start agagin as the threshold increased
                }
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

        size_t corrected_estimate_per_part = estimate_per_part * 1.5;// * (1.0 + 2 * static_cast<double>(joint_estimate)/sum_of_cardinalities);

        size_t const u_bins_per_part = seqan::hibf::divide_and_ceil(positions.size(), config.number_of_partitions);
        size_t const block_size =
            std::min(u_bins_per_part,
                     chopper::next_multiple_of_64(static_cast<uint16_t>(std::ceil(std::sqrt(u_bins_per_part)))));
        size_t const number_of_blocks = seqan::hibf::divide_and_ceil(positions.size(), block_size);

        std::vector<size_t> max_partition_cardinality(config.number_of_partitions, 0u);
        std::vector<size_t> min_partition_cardinality(config.number_of_partitions, std::numeric_limits<size_t>::max());

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
                size_t const user_bin_idx{sorted_positions[block_size * i + x]};
                partition_sketches[p].merge(sketches[user_bin_idx]);
                partitions[p].push_back(user_bin_idx);
                max_partition_cardinality[p] = std::max(max_partition_cardinality[p], cardinalities[user_bin_idx]);
                min_partition_cardinality[p] = std::min(min_partition_cardinality[p], cardinalities[user_bin_idx]);
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
                                corrected_estimate_per_part,
                                {sorted_positions[i]},
                                cardinalities,
                                sketches,
                                partitions,
                                partition_sketches,
                                max_partition_cardinality,
                                min_partition_cardinality);
        }
    }
    else if (config.partitioning_approach == partitioning_scheme::lsh)
    {
        std::vector<size_t> partition_sizes(config.number_of_partitions, 0u);

        size_t corrected_estimate_per_part = seqan::hibf::divide_and_ceil(sum_of_cardinalities, config.number_of_partitions); //estimate_per_part * 1.5;// * (1.0 + 2 * static_cast<double>(joint_estimate)/sum_of_cardinalities);
        size_t original_estimate_per_part = estimate_per_part;

        std::vector<size_t> max_partition_cardinality(config.number_of_partitions, 0u);
        std::vector<size_t> min_partition_cardinality(config.number_of_partitions, std::numeric_limits<size_t>::max());

        // initial partitioning using locality sensitive hashing (LSH)
        config.lsh_algorithm_timer.start();
        std::vector<Cluster> clusters = very_similar_LSH_partitioning(minHash_sketches, positions, cardinalities, sketches, original_estimate_per_part, config);
        post_process_clusters(clusters, cardinalities, config);
        config.lsh_algorithm_timer.stop();

std::ofstream ofs{"/tmp/final.clusters"};
for (size_t i = 0; i < clusters.size(); ++i)
{
    seqan::hibf::sketch::hyperloglog sketch(config.hibf_config.sketch_bits);
    for (size_t j = 0; j < clusters[i].size(); ++j)
    {
        sketch.merge(sketches[clusters[i].contained_user_bins()[j]]);
    }
    ofs << i << ":" << sketch.estimate() << ":" << clusters[i].size() << std::endl;
}

        // initialise partitions with the first p largest clusters (post_processing sorts by size)
        size_t cidx{0}; // current cluster index
        for (size_t p = 0; p < config.number_of_partitions; ++p)
        {
            assert(!clusters[cidx].empty());
            auto const & cluster = clusters[cidx].contained_user_bins();
            bool split_cluster = false;

            if (cluster.size() > config.hibf_config.tmax)
            {
                size_t card{0};
                for (size_t uidx = 0; uidx < cluster.size(); ++uidx)
                    card += cardinalities[cluster[uidx]];

                if (card > 0.05 * sum_of_cardinalities / config.hibf_config.tmax)
                    split_cluster = true;
            }

            size_t const end = (split_cluster) ? cluster.size() : std::min(cluster.size(), config.hibf_config.tmax);
            for (size_t uidx = 0; uidx < end; ++uidx)
            {
                size_t const user_bin_idx = cluster[uidx];
                // if a single cluster already exceeds the cardinality_per_part,
                // then the remaining user bins of the cluster must spill over into the next partition
                if ((uidx != 0 && (uidx % config.hibf_config.tmax == 0)) || partition_sizes[p] > corrected_estimate_per_part)
                {
                    ++p;
                    // p_has_been_incremented = true;
                    assert(p < config.number_of_partitions);
                }

                partition_sizes[p] += cardinalities[user_bin_idx];
                partitions[p].push_back(user_bin_idx);
                max_partition_cardinality[p] = std::max(max_partition_cardinality[p], cardinalities[user_bin_idx]);
                min_partition_cardinality[p] = std::min(min_partition_cardinality[p], cardinalities[user_bin_idx]);
            }

            ++cidx;
        }

        std::vector<std::vector<size_t>> remaining_clusters{};

        for (size_t i = 0; i < clusters.size(); ++i)
        {
            if (clusters[i].empty())
                break;

            auto const & cluster = clusters[i].contained_user_bins();

            // if i < cidx the first cluster was already used in initialising above
            if (i < cidx)
            {
                bool split_cluster = false;
                if (cluster.size() > config.hibf_config.tmax)
                {
                    size_t card{0};
                    for (size_t uidx = 0; uidx < cluster.size(); ++uidx)
                        card += cardinalities[cluster[uidx]];

                    if (card > 0.05 * sum_of_cardinalities / config.hibf_config.tmax)
                        split_cluster = true;
                }

                if (cluster.size() > config.hibf_config.tmax && !split_cluster) // then only the first tmax UBs have been consumed
                {
                    std::vector<size_t> remainder(cluster.begin() + config.hibf_config.tmax, cluster.end());
                    remaining_clusters.insert(remaining_clusters.end(), remainder);
                }
                // otherwise ignore cluster
            }
            else
            {
                remaining_clusters.insert(remaining_clusters.end(), cluster);
            }
        }

        // assign the rest by similarity
        for (size_t ridx = 0; ridx < remaining_clusters.size(); ++ridx)
        {
            auto const & cluster = remaining_clusters[ridx];

            config.search_partition_algorithm_timer.start();
            find_best_sorted_partition(config, corrected_estimate_per_part, cluster, cardinalities, sketches, partitions, partition_sizes, max_partition_cardinality, min_partition_cardinality);
            config.search_partition_algorithm_timer.start();

        }
    }
    else if (config.partitioning_approach == partitioning_scheme::lsh_sim)
    {
        uint8_t const sketch_bits{config.hibf_config.sketch_bits};
//std::cout << "LSH partitioning into " << config.number_of_partitions << std::endl;
        std::vector<seqan::hibf::sketch::hyperloglog> partition_sketches(config.number_of_partitions,
                                                                         seqan::hibf::sketch::hyperloglog(sketch_bits));

        size_t corrected_estimate_per_part = estimate_per_part * 1.5;// * (1.0 + 2 * static_cast<double>(joint_estimate)/sum_of_cardinalities);
        size_t original_estimate_per_part = estimate_per_part;

        std::vector<size_t> max_partition_cardinality(config.number_of_partitions, 0u);
        std::vector<size_t> min_partition_cardinality(config.number_of_partitions, std::numeric_limits<size_t>::max());

        // initial partitioning using locality sensitive hashing (LSH)
        config.lsh_algorithm_timer.start();
        std::vector<Cluster> clusters = very_similar_LSH_partitioning(minHash_sketches, positions, cardinalities, sketches, original_estimate_per_part, config);
        post_process_clusters(clusters, cardinalities, config);
        config.lsh_algorithm_timer.stop();

std::ofstream ofs{"/tmp/final.clusters"};
for (size_t i = 0; i < clusters.size(); ++i)
{
    seqan::hibf::sketch::hyperloglog sketch(config.hibf_config.sketch_bits);
    for (size_t j = 0; j < clusters[i].size(); ++j)
    {
        sketch.merge(sketches[clusters[i].contained_user_bins()[j]]);
    }
    ofs << i << ":" << sketch.estimate() << ":" << clusters[i].size() << std::endl;
}

        // initialise partitions with the first p largest clusters (post_processing sorts by size)
        size_t cidx{0}; // current cluster index
        for (size_t p = 0; p < config.number_of_partitions; ++p)
        {
            assert(!clusters[cidx].empty());
            auto const & cluster = clusters[cidx].contained_user_bins();
            bool split_cluster = false;

            if (cluster.size() > config.hibf_config.tmax)
            {
                size_t card{0};
                for (size_t uidx = 0; uidx < cluster.size(); ++uidx)
                    card += cardinalities[cluster[uidx]];

                if (card > 0.05 * sum_of_cardinalities / config.hibf_config.tmax)
                    split_cluster = true;
            }

            size_t const end = (split_cluster) ? cluster.size() : std::min(cluster.size(), config.hibf_config.tmax);
            for (size_t uidx = 0; uidx < end; ++uidx)
            {
                size_t const user_bin_idx = cluster[uidx];
                // if a single cluster already exceeds the cardinality_per_part,
                // then the remaining user bins of the cluster must spill over into the next partition
                if ((uidx != 0 && (uidx % config.hibf_config.tmax == 0)) || partition_sketches[p].estimate() > corrected_estimate_per_part)
                {
                    ++p;
                    // p_has_been_incremented = true;
                    assert(p < config.number_of_partitions);
                }

                partition_sketches[p].merge(sketches[user_bin_idx]);
                partitions[p].push_back(user_bin_idx);
                max_partition_cardinality[p] = std::max(max_partition_cardinality[p], cardinalities[user_bin_idx]);
                min_partition_cardinality[p] = std::min(min_partition_cardinality[p], cardinalities[user_bin_idx]);
            }

            ++cidx;
        }

        std::vector<std::vector<size_t>> remaining_clusters{};

        for (size_t i = 0; i < clusters.size(); ++i)
        {
            if (clusters[i].empty())
                break;

            auto const & cluster = clusters[i].contained_user_bins();

            // if i < cidx the first cluster was already used in initialising above
            if (i < cidx)
            {
                bool split_cluster = false;
                if (cluster.size() > config.hibf_config.tmax)
                {
                    size_t card{0};
                    for (size_t uidx = 0; uidx < cluster.size(); ++uidx)
                        card += cardinalities[cluster[uidx]];

                    if (card > 0.05 * sum_of_cardinalities / config.hibf_config.tmax)
                        split_cluster = true;
                }

                if (cluster.size() > config.hibf_config.tmax && !split_cluster) // then only the first tmax UBs have been consumed
                {
                    std::vector<size_t> remainder(cluster.begin() + config.hibf_config.tmax, cluster.end());
                    remaining_clusters.insert(remaining_clusters.end(), remainder);
                }
                // otherwise ignore cluster
            }
            else
            {
                remaining_clusters.insert(remaining_clusters.end(), cluster);
            }
        }

        // assign the rest by similarity
        for (size_t ridx = 0; ridx < remaining_clusters.size(); ++ridx)
        {
            auto const & cluster = remaining_clusters[ridx];

            config.search_partition_algorithm_timer.start();
            find_best_partition(config, corrected_estimate_per_part, cluster, cardinalities, sketches, partitions, partition_sketches, max_partition_cardinality, min_partition_cardinality);
            config.search_partition_algorithm_timer.start();
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

    // assert([&](){ bool x{false}; for (auto const & p : partitions) { x &= !p.empty(); }; return x; });
}

seqan::hibf::layout::layout general_layout(chopper::configuration const & config,
                std::vector<size_t> positions,
                std::vector<size_t> const & cardinalities,
                std::vector<seqan::hibf::sketch::hyperloglog> const & sketches)
{
    seqan::hibf::layout::layout hibf_layout;

    seqan::hibf::concurrent_timer union_estimation_timer{};
    seqan::hibf::concurrent_timer rearrangement_timer{};
    seqan::hibf::concurrent_timer dp_algorithm_timer{};

    dp_algorithm_timer.start();
    hibf_layout = seqan::hibf::layout::compute_layout(config.hibf_config,
                                                      cardinalities,
                                                      sketches,
                                                      std::move(positions),
                                                      union_estimation_timer,
                                                      rearrangement_timer);
    dp_algorithm_timer.stop();

    return hibf_layout;
}

bool do_I_need_a_fast_layout(chopper::configuration const & config, std::vector<size_t> const & positions, std::vector<size_t> const & cardinalities)
{
    // the fast layout heuristic would greedily merge even if merging only 2 bins at a time
    // merging only little number of bins is highly disadvantegous for lower levels because few bins
    // will be heavily split and this will raise the fpr correction for split bins
    // Thus, if the average number of user bins per technical bin is less then 64, we should not fast layout
    if (positions.size() < (64 * config.hibf_config.tmax))
        return false;

    if (positions.size() > 500'000) // layout takes more than half a day (should this be a user option?)
        return true;

    size_t largest_size{0};
    size_t sum_of_cardinalities{0};

    for (size_t const i : positions)
    {
        sum_of_cardinalities += cardinalities[i];
        largest_size = std::max(largest_size, cardinalities[i]);
    }

    size_t const cardinality_per_tb = sum_of_cardinalities / config.hibf_config.tmax;

    bool const largest_user_bin_might_be_split = largest_size > cardinality_per_tb;

    // if no splitting is needed, its worth it to use a fast-merge-only algorithm
    if (!largest_user_bin_might_be_split)
        return true;

    return false;
}

size_t determine_max_bin(std::vector<std::vector<size_t>> const & positions,
                         std::vector<seqan::hibf::sketch::hyperloglog> const & sketches)
{
    size_t max_bin_id{0};
    size_t max_size{0};

    for (size_t i = 0; i < positions.size(); ++i)
    {
        if (positions[i].empty()) // should only happen on really small data sets with edge cases
            continue;

        seqan::hibf::sketch::hyperloglog tb{sketches[positions[i][0]]}; // init
        for (auto const pos : positions[i])
            tb.merge(sketches[pos]);

        if (tb.estimate() > max_size)
        {
            max_bin_id = i;
            max_size = tb.estimate();
        }
    }

    return max_bin_id;
}

void add_level_to_layout(seqan::hibf::layout::layout & hibf_layout,
                         std::vector<std::vector<size_t>> const & positions,
                         std::vector<seqan::hibf::sketch::hyperloglog> const & sketches,
                         std::vector<size_t> const & previous)
{
    hibf_layout.max_bins.emplace_back(previous, determine_max_bin(positions, sketches)); // add lower level meta information

    // we assume here that the user bins have been sorted by user bin id such that pos = idx
    size_t partition_idx{0};
    for (auto const & partition : positions)
    {
        for (size_t const user_bin_id : partition)
        {
            assert(hibf_layout.user_bins[user_bin_id].idx == user_bin_id);
            auto & current_user_bin = hibf_layout.user_bins[user_bin_id];

            // update
            assert(previous == current_user_bin.previous_TB_indices);
            // TODO: this should not happen for single bins?
            current_user_bin.previous_TB_indices.push_back(partition_idx);
        }
        ++partition_idx;
    }
}

void update_layout_from_child_layout(seqan::hibf::layout::layout & child_layout,
                                     seqan::hibf::layout::layout & hibf_layout,
                                     std::vector<size_t> const & new_previous)
{
    hibf_layout.max_bins.emplace_back(new_previous, child_layout.top_level_max_bin_id);

    for (auto & max_bin : child_layout.max_bins)
    {
        max_bin.previous_TB_indices.insert(max_bin.previous_TB_indices.begin(), new_previous.begin(), new_previous.end());
        hibf_layout.max_bins.push_back(max_bin);
    }

    for (auto const & user_bin : child_layout.user_bins)
    {
        auto & actual_user_bin = hibf_layout.user_bins[user_bin.idx];

        actual_user_bin.previous_TB_indices.insert(actual_user_bin.previous_TB_indices.end(),
                                                    user_bin.previous_TB_indices.begin(),
                                                    user_bin.previous_TB_indices.end());
        actual_user_bin.number_of_technical_bins = user_bin.number_of_technical_bins;
        actual_user_bin.storage_TB_id = user_bin.storage_TB_id;
    }
}

void fast_layout_recursion(chopper::configuration const & config,
            std::vector<size_t> const & positions,
            std::vector<size_t> const & cardinalities,
            std::vector<seqan::hibf::sketch::hyperloglog> const & sketches,
            std::vector<seqan::hibf::sketch::minhashes> const & minHash_sketches,
            seqan::hibf::layout::layout & hibf_layout,
            std::vector<size_t> const & previous)
{
    std::vector<std::vector<size_t>> tmax_partitions(config.hibf_config.tmax);

    // here we assume that we want to start with a fast layout
    partition_user_bins(config, positions, cardinalities, sketches, minHash_sketches, tmax_partitions);

    #pragma omp critical
    {
        add_level_to_layout(hibf_layout, tmax_partitions, sketches, previous);
    }

    for (size_t partition_idx = 0; partition_idx < tmax_partitions.size(); ++partition_idx)
    {
        auto const & partition = tmax_partitions[partition_idx];
        auto const new_previous = [&] () {auto cpy{previous}; cpy.push_back(partition_idx); return cpy; }();

        if (partition.empty() && partition.size() == 1) // nothing to merge
            continue;

        if (do_I_need_a_fast_layout(config, partition, cardinalities))
        {
            fast_layout_recursion(config, partition, cardinalities, sketches, minHash_sketches, hibf_layout, previous); // recurse fast_layout
        }
        else
        {
            auto child_layout = general_layout(config, partition, cardinalities, sketches);

            #pragma omp critical
            {
                update_layout_from_child_layout(child_layout, hibf_layout, new_previous);
            }
        }
    }
}

void fast_layout(chopper::configuration const & config,
                std::vector<size_t> const & positions,
                std::vector<size_t> const & cardinalities,
                std::vector<seqan::hibf::sketch::hyperloglog> const & sketches,
                std::vector<seqan::hibf::sketch::minhashes> const & minHash_sketches,
                seqan::hibf::layout::layout & hibf_layout)
{
    auto config_copy = config;
    config_copy.number_of_partitions = config.hibf_config.tmax;

    std::vector<std::vector<size_t>> tmax_partitions(config.hibf_config.tmax);

    // here we assume that we want to start with a fast layout
    config.intital_partition_timer.start();
    partition_user_bins(config_copy, positions, cardinalities, sketches, minHash_sketches, tmax_partitions);
    config.intital_partition_timer.stop();

    hibf_layout.top_level_max_bin_id = determine_max_bin(tmax_partitions, sketches);

    hibf_layout.user_bins.resize(config_copy.hibf_config.number_of_user_bins);

    // initialise user bins in layout
    for (size_t partition_idx = 0; partition_idx < tmax_partitions.size(); ++partition_idx)
    {
        for (size_t const user_bin_id : tmax_partitions[partition_idx])
        {
            hibf_layout.user_bins[user_bin_id] = {.previous_TB_indices = {partition_idx},
                                                  .storage_TB_id = 0 /*not determiend yet*/,
                                                  .number_of_technical_bins = 1 /*not determiend yet*/,
                                                  .idx = user_bin_id};
        }
    }

    config.small_layouts_timer.start();
    #pragma omp parallel
    #pragma omp single
    {
        #pragma omp taskloop
        for (size_t partition_idx = 0; partition_idx < tmax_partitions.size(); ++partition_idx)
        {
            auto const & partition = tmax_partitions[partition_idx];

            if (partition.empty() || partition.size() == 1) // nothing to merge
                continue;

            if (do_I_need_a_fast_layout(config_copy, partition, cardinalities))
            {
                fast_layout_recursion(config_copy, partition, cardinalities, sketches, minHash_sketches, hibf_layout, {partition_idx}); // recurse fast_layout
            }
            else
            {
                auto small_layout = general_layout(config, partition, cardinalities, sketches);

                #pragma omp critical
                {
                    update_layout_from_child_layout(small_layout, hibf_layout, std::vector<size_t>{partition_idx});
                }
            }
        }
    }
    config.small_layouts_timer.stop();
}

int execute(chopper::configuration & config,
            std::vector<std::vector<std::string>> const & filenames,
            std::vector<seqan::hibf::sketch::hyperloglog> const & sketches,
            std::vector<seqan::hibf::sketch::minhashes> const & minHash_sketches)
{
    config.hibf_config.validate_and_set_defaults();

    std::vector<size_t> cardinalities;
    seqan::hibf::sketch::estimate_kmer_counts(sketches, cardinalities);

    if (config.number_of_partitions < 2) // 0 == unset == single HIBF, 1 == single HIBF
    {
        seqan::hibf::layout::layout hibf_layout;


        if (config.determine_best_tmax)
        {
            hibf_layout = determine_best_number_of_technical_bins(config, cardinalities, sketches);
        }
        else
        {
            config.dp_algorithm_timer.start();
            fast_layout(config, seqan::hibf::iota_vector(sketches.size()), cardinalities, sketches, minHash_sketches, hibf_layout);
            config.dp_algorithm_timer.stop();

            // hibf_layout = seqan::hibf::layout::compute_layout(config.hibf_config,
            //                                                   cardinalities,
            //                                                   sketches,
            //                                                   seqan::hibf::iota_vector(sketches.size()),
            //                                                   config.union_estimation_timer,
            //                                                   config.rearrangement_timer);

            // sort records ascending by the number of bin indices (corresponds to the IBF levels)
            // GCOVR_EXCL_START
            std::ranges::sort(hibf_layout.max_bins,
                            [](auto const & r, auto const & l)
                            {
                                return r.previous_TB_indices.size() < l.previous_TB_indices.size();
                            });
            // GCOVR_EXCL_STOP

            if (config.output_verbose_statistics)
            {
                size_t dummy{};
                chopper::layout::hibf_statistics global_stats{config, sketches, cardinalities};
                global_stats.hibf_layout = hibf_layout;
                global_stats.print_header_to(std::cout);
                global_stats.print_summary_to(dummy, std::cout);
            }
        }

        // brief Write the output to the layout file.
        std::ofstream fout{config.output_filename};
        chopper::layout::write_user_bins_to(filenames, fout);
        config.write_to(fout);
        hibf_layout.write_to(fout);
    }
    else
    {
        std::vector<std::vector<size_t>> positions(config.number_of_partitions); // asign positions for each partition

        std::vector<size_t> ps;
        ps.resize(cardinalities.size());
        std::iota(ps.begin(), ps.end(), 0);
        partition_user_bins(config, ps, cardinalities, sketches, minHash_sketches, positions);

        std::vector<seqan::hibf::layout::layout> hibf_layouts(config.number_of_partitions); // multiple layouts

#pragma omp parallel for schedule(dynamic) num_threads(config.hibf_config.threads)
        for (size_t i = 0; i < config.number_of_partitions; ++i)
        {
            // reset tmax to fit number of user bins in layout
            auto local_hibf_config = config.hibf_config; // every thread needs to set individual tmax
            local_hibf_config.tmax =
                chopper::next_multiple_of_64(static_cast<uint16_t>(std::ceil(std::sqrt(positions[i].size()))));

            config.dp_algorithm_timer.start();
            hibf_layouts[i] = seqan::hibf::layout::compute_layout(local_hibf_config,
                                                                  cardinalities,
                                                                  sketches,
                                                                  std::move(positions[i]),
                                                                  config.union_estimation_timer,
                                                                  config.rearrangement_timer);
            config.dp_algorithm_timer.stop();
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
