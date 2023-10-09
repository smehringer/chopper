// ---------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
// ---------------------------------------------------------------------------------------------------

#include <iostream>
#include <set>

#include <robin_hood.h>

#include <sharg/detail/to_string.hpp>
#include <sharg/exceptions.hpp>
#include <sharg/parser.hpp>

#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/search/views/minimiser_hash.hpp>

#include <chopper/adjust_seed.hpp>
#include <chopper/layout/hibf_statistics.hpp>
#include <chopper/layout/input.hpp>

#include <hibf/build/bin_size_in_bits.hpp>
#include <hibf/build/build_data.hpp>
#include <hibf/interleaved_bloom_filter.hpp>
#include <hibf/layout/compute_fpr_correction.hpp>
#include <hibf/layout/graph.hpp>
#include <hibf/sketch/hyperloglog.hpp>

#include "shared.hpp"

void init_shared_meta(sharg::parser & parser)
{
    parser.info.version = "1.0.0";
    parser.info.author = "Svenja Mehringer";
    parser.info.email = "svenja.mehringer@fu-berlin.de";
    parser.info.short_description = "Compute a top-level HIBF layout figure file";
    parser.info.description.emplace_back("Computes an table to display the top-level layout.");
}

void init_options(sharg::parser & parser, config & cfg)
{
    parser.add_subsection("Main options:");
    parser.add_option(
        cfg.input,
        sharg::config{.short_id = '\0',
                      .long_id = "input",
                      .description = "The input must be a layout file computed via chopper layout or raptor layout. ",
                      .required = true,
                      .validator = sharg::input_file_validator{}});
    parser.add_option(cfg.output,
                      sharg::config{.short_id = '\0',
                                    .long_id = "output",
                                    .description = "The output. ",
                                    .required = true,
                                    .validator = sharg::output_file_validator{}});
    parser.add_option(
        cfg.threads,
        sharg::config{.short_id = '\0', .long_id = "threads", .description = "The number of threads to use."});
}

void parse(sharg::parser & parser)
{
    try
    {
        parser.parse();
    }
    catch (sharg::parser_error const & ext)
    {
        std::cerr << "[ERROR] " << ext.what() << '\n';
        std::exit(EXIT_FAILURE);
    }
}

int main(int argc, char const * argv[])
{
    config cfg{};

    sharg::parser main_parser{"layout_stats", argc, argv, sharg::update_notifications::off, {"general", "sizes"}};
    init_shared_meta(main_parser);

    parse(main_parser);
    sharg::parser & sub_parser = main_parser.get_sub_parser();

    init_shared_meta(sub_parser);
    init_options(sub_parser, cfg);
    parse(sub_parser);

    if (sub_parser.info.app_name == std::string_view{"layout_stats-general"})
        execute_general(cfg);
    else if (sub_parser.info.app_name == std::string_view{"layout_stats-sizes"})
        execute_sizes(cfg);
    else
        std::cerr << "[ERROR] Unknown subcommand\n";
}
