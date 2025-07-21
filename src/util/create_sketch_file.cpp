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
    sharg::parser parser{"create_sketch_file", argc, argv, sharg::update_notifications::off};
    parser.info.version = "1.0.0";

    sketch_file sout{};

    std::filesystem::path output_file{"out.sketches"};
    set_up_parser(parser, sout.chopper_config, output_file);

    try
    {
        parser.parse();
    }
    catch (sharg::parser_error const & ext)
    {
        std::cerr << "[ERROR] " << ext.what() << '\n';
        return -1;
    }

    chopper::sketch::read_data_file(sout.chopper_config, sout.filenames);

    sout.chopper_config.hibf_config.input_fn =
        chopper::input_functor{sout.filenames, sout.chopper_config.precomputed_files, sout.chopper_config.k, sout.chopper_config.window_size};
    sout.chopper_config.hibf_config.number_of_user_bins = sout.filenames.size();

    seqan::hibf::sketch::compute_sketches(sout.chopper_config.hibf_config, sout.hll_sketches, sout.minHash_sketches);

    std::ofstream os{output_file, std::ios::binary};
    cereal::BinaryOutputArchive oarchive{os};
    oarchive(sout);
}
