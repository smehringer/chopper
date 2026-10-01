// ---------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/chopper/blob/main/LICENSE.md
// ---------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides chopper::layout::phibf::execute.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <cstddef>
#include <string>
#include <vector>

#include <chopper/configuration.hpp>

#include <hibf/sketch/hyperloglog.hpp>
#include <hibf/sketch/minhashes.hpp>

namespace chopper::layout::phibf
{

/*!\brief Computes a partitioned HIBF layout and writes it to `config.output_filename`.
 * \param[in,out] config           The configuration. `config.number_of_partitions` must be at least 2.
 *                                 `config.hibf_config` must be validated (`validate_and_set_defaults`). The timers are
 *                                 updated.
 * \param[in]     filenames        The file names of each user bin. They are written to the layout file.
 * \param[in]     cardinalities    The cardinality of each user bin.
 * \param[in]     sketches         The HyperLogLog sketch of each user bin.
 * \param[in]     minHash_sketches The MinHash sketches of each user bin. Required by the `lsh` and `lsh_sim`
 *                                 partitioning approaches.
 * \returns 0.
 *
 * 1. The user bins are distributed onto `config.number_of_partitions` partitions (phibf::partition_user_bins).
 * 2. Each partition is laid out (in parallel) with the DP layout of the HIBF library
 *    (seqan::hibf::layout::compute_layout), with `tmax` set to `next_multiple_of_64(ceil(sqrt(#user bins)))` of the
 *    partition.
 * 3. The layout file contains the user bins and the configuration once, followed by one layout per partition.
 *    The user bin indices in each layout are global. Read it with chopper::layout::read_layouts_file.
 */
int execute(chopper::configuration & config,
            std::vector<std::vector<std::string>> const & filenames,
            std::vector<size_t> const & cardinalities,
            std::vector<seqan::hibf::sketch::hyperloglog> const & sketches,
            std::vector<seqan::hibf::sketch::minhashes> const & minHash_sketches);

} // namespace chopper::layout::phibf
