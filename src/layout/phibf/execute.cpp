// ---------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
// ---------------------------------------------------------------------------------------------------

#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <string>
#include <utility>
#include <vector>

#include <chopper/configuration.hpp>
#include <chopper/layout/output.hpp>
#include <chopper/layout/phibf/execute.hpp>
#include <chopper/layout/phibf/partition_user_bins.hpp>
#include <chopper/next_multiple_of_64.hpp>

#include <hibf/layout/compute_layout.hpp>
#include <hibf/layout/layout.hpp>

namespace chopper::layout::phibf
{

int execute(chopper::configuration & config,
            std::vector<std::vector<std::string>> const & filenames,
            std::vector<size_t> const & cardinalities,
            std::vector<seqan::hibf::sketch::hyperloglog> const & sketches,
            std::vector<seqan::hibf::sketch::minhashes> const & minHash_sketches)
{
    assert(config.number_of_partitions >= 2u);

    std::vector<std::vector<size_t>> positions(config.number_of_partitions); // asign positions for each partition

    partition_user_bins(config, cardinalities, sketches, minHash_sketches, positions);

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

    return 0;
}

} // namespace chopper::layout::phibf
