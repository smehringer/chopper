#pragma once

#include <seqan3/std/filesystem>
#include <fstream>
#include <future>
#include <thread>

#include <robin_hood.h>

#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/io/views/async_input_buffer.hpp>
#include <seqan3/search/views/kmer_hash.hpp>
#include <seqan3/search/views/minimiser_hash.hpp>
#include <seqan3/utility/views/to.hpp>

#include <chopper/count/configuration.hpp>
#include <chopper/count/output.hpp>
#include <chopper/sketch/hyperloglog.hpp>

namespace chopper::count
{

struct mytraits : public seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = seqan3::dna4;
};

using sequence_file_type = seqan3::sequence_file_input<mytraits,
                                                       seqan3::fields<seqan3::field::seq>,
                                                       seqan3::type_list<seqan3::format_fasta, seqan3::format_fastq>>;

template <typename hash_view_type>
inline void count_kmers(robin_hood::unordered_map<std::string, std::vector<std::string>> const & filename_clusters,
                        configuration const & config,
                        hash_view_type && hash_fn)
{
   // output file
    std::ofstream fout{config.count_filename};

    if (!fout.good())
        throw std::runtime_error{"Could not open file" + config.count_filename.string() + " for reading."};

    // create the hll dir if it doesn't already exist
    if (!config.disable_sketch_output)
        std::filesystem::create_directory(config.sketch_directory);

    // copy filename clusters to vector
    std::vector<std::pair<std::string, std::vector<std::string>>> cluster_vector{};
    for (auto const & cluster : filename_clusters)
        cluster_vector.emplace_back(cluster.first, cluster.second);

    #pragma omp parallel for schedule(static) num_threads(config.num_threads)
    for (size_t i = 0; i < cluster_vector.size(); ++i)
    {
        chopper::sketch::hyperloglog sketch(config.sketch_bits);

        // read files
        for (auto const & filename : cluster_vector[i].second)
            for (auto && [seq] : sequence_file_type{filename})
                for (auto && hash : seq | hash_fn)
                    sketch.add(reinterpret_cast<char*>(&hash), sizeof(hash));

        // print either the exact or the approximate count, depending on exclusively_hlls
        uint64_t const weight = sketch.estimate();

        #pragma omp critical
        write_count_file_line(cluster_vector[i], weight, fout);

        if (!config.disable_sketch_output)
            write_sketch_file(cluster_vector[i], sketch, config);
    }
}

inline void count_kmers(robin_hood::unordered_map<std::string, std::vector<std::string>> const & filename_clusters,
                        configuration const & config)
{

    auto compute_minimiser = seqan3::views::minimiser_hash(seqan3::ungapped{config.k},
                                                           seqan3::window_size{config.w},
                                                           seqan3::seed{0x8F3F73B5CF1C9ADE >> (64u - 2u * config.k)});

    auto compute_kmers = seqan3::views::kmer_hash(seqan3::ungapped{config.k});

    if (config.disable_minimizers)
        count_kmers(filename_clusters, config, compute_kmers);
    else
        count_kmers(filename_clusters, config, compute_minimiser);
}

} // namespace chopper::count
