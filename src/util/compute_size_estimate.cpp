// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <iostream>

#include <sharg/detail/to_string.hpp>
#include <sharg/exceptions.hpp>
#include <sharg/parser.hpp>

#include <chopper/configuration.hpp>
#include <chopper/input_functor.hpp>
#include <chopper/sketch/sketch_file.hpp>
#include <chopper/sketch/read_data_file.hpp>

#include <hibf/sketch/compute_sketches.hpp>
#include <hibf/sketch/hyperloglog.hpp>
#include <hibf/sketch/estimate_kmer_counts.hpp> // for estimate_kmer_counts

inline void set_up_parser(sharg::parser & parser, chopper::configuration & config, std::filesystem::path & out_file)
{
    parser.info.version = "1.0.0";
    parser.info.author = "Svenja Mehringer";
    parser.info.email = "svenja.mehringer@fu-berlin.de";
    parser.info.short_description = "Read IBF sizes from an HIBF";

    parser.add_subsection("Main options:");
    parser.add_option(
        config.data_file,
        sharg::config{
            .short_id = '\0',
            .long_id = "input",
            .description =
                "The input must be a file containing paths to sequence data you wish to estimate; one filepath "
                "per line. If your file contains auxiliary information (e.g. species IDs), your file must be tab-"
                "separated.",
            .required = true});
    parser.add_list_item("", "Example file:");
    parser.add_list_item("", "```");
    parser.add_list_item("", "/absolute/path/to/file1.fasta");
    parser.add_list_item("", "/absolute/path/to/file2.fa.gz");
    parser.add_list_item("", "```");

    parser.add_option(out_file,
                      sharg::config{.short_id = '\0',
                                    .long_id = "output",
                                    .description = "The output. ",
                                    .required = true,
                                    .validator = sharg::output_file_validator{}});

    parser.add_option(
        config.k,
        sharg::config{
            .short_id = '\0',
            .long_id = "kmer",
            .description =
                "The k-mer size influences the size estimates of the input. "
                "Choosing a k-mer size that is too small for "
                "your data will result in files appearing more similar than they really are. Likewise, a large "
                "k-mer size might miss out on certain similarities. For DNA sequences, a k-mer size between "
                "[16,32] has proven to work well."});

    parser.add_option(
        config.window_size,
        sharg::config{
            .short_id = '\0',
            .long_id = "window",
            .description =
                "Setting this option will trigger the computation of (w,k)-minimizers instead of canonical kmers. "
                "Minimizers can thin out the data, reduce the memory footpring of the resulting index and increase "
                "runtime performance. On the other hand, it also decreases accuracy and might cause false negatives."
                "For DNA sequences, a window size of 2-4 positions more than the kmers size, e.g. (42,20)-minimizers, "
                "has proven to reduce the computational effort significantly while only slightly descreasing the "
                "accuracy.",
            .default_message = "k-mer size",
        });

    parser.add_option(
        config.hibf_config.sketch_bits,
        sharg::config{.short_id = '\0',
                      .long_id = "sketch-bits",
                      .description =
                          "The number of bits the HyperLogLog sketch should use to distribute the values into bins.",
                      .advanced = true,
                      .validator = sharg::arithmetic_range_validator{5, 32}});

    parser.add_option(
        config.hibf_config.threads,
        sharg::config{
            .short_id = '\0',
            .long_id = "threads",
            .description =
                "The number of threads to use. Currently, only merging of sketches is parallelized, so if the flag "
                "--disable-rearrangement is set, --threads will have no effect.",
            .validator =
                sharg::arithmetic_range_validator{static_cast<size_t>(1), std::numeric_limits<size_t>::max()}});
}

int main(int argc, char const * argv[])
{
    sharg::parser parser{"compute_estimates", argc, argv, sharg::update_notifications::off};
    parser.info.version = "1.0.0";

    chopper::configuration chopper_config{};
    std::vector<std::vector<std::string>> filenames{};
    std::vector<seqan::hibf::sketch::hyperloglog> sketches{};
    std::vector<seqan::hibf::sketch::minhashes> minHash_sketches{};
    std::vector<size_t> cardinalities;

    std::filesystem::path output_file{"out.estimates"};
    set_up_parser(parser, chopper_config, output_file);

    try
    {
        parser.parse();
    }
    catch (sharg::parser_error const & ext)
    {
        std::cerr << "[ERROR] " << ext.what() << '\n';
        return -1;
    }

    auto has_sketch_file_extension = [](std::filesystem::path const & path)
    {
        return path.string().ends_with(".sketch") || path.string().ends_with(".sketches");
    };

    bool const input_is_a_sketch_file = has_sketch_file_extension(chopper_config.data_file);

    if (input_is_a_sketch_file)
    {
        chopper::sketch::sketch_file sin{};

        { // Deserialization is guaranteed to be complete when going out of scope.
            std::ifstream is{chopper_config.data_file};
            cereal::BinaryInputArchive iarchive{is};
            iarchive(sin);
        }

        filenames = std::move(sin.filenames); // No need to call check_filenames because the files are not read.
        sketches = std::move(sin.hll_sketches);
    }
    else
    {
        chopper::sketch::read_data_file(chopper_config, filenames);

        chopper_config.hibf_config.input_fn =
            chopper::input_functor{filenames, chopper_config.precomputed_files, chopper_config.k, chopper_config.window_size};
        chopper_config.hibf_config.number_of_user_bins = filenames.size();

        seqan::hibf::sketch::compute_sketches(chopper_config.hibf_config, sketches, minHash_sketches);
    }

    seqan::hibf::sketch::estimate_kmer_counts(sketches, cardinalities);

    // compute union estimate
    auto union_sketch = sketches[0];
    for (size_t i{0}; i < sketches.size(); ++i)
        union_sketch.merge(sketches[i]);

    std::cout << "union estimate of all sketches merged: " << union_sketch.estimate() << std::endl;

    std::ofstream fout{output_file};

    for (size_t i{0}; i < filenames.size(); ++i)
    {
        fout << filenames[i][0] << "\t" << cardinalities[i] << std::endl;
    }
}
