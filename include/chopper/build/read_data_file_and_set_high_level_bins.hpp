#pragma once

#include <set>
#include <string>
#include <vector>

#include <seqan3/std/charconv>

#include <chopper/build/build_config.hpp>
#include <chopper/build/header_data.hpp>
#include <chopper/detail_bin_prefixes.hpp>
#include <chopper/detail_parse_binning_line.hpp>
#include <chopper/detail_starts_with.hpp>

auto read_data_file_and_set_high_level_bins(build_config const & config)
{
    header_data header;
    std::vector<data_file_record> records{};

    std::unordered_map<std::string, data_file_record> low_level_records{};

    std::ifstream binning_file{config.binning_filename};

    if (!binning_file.good() || !binning_file.is_open())
        throw std::logic_error{"Could not open file for reading"};

    std::string current_line;
    while (std::getline(binning_file, current_line) && current_line[0] == '#')
    {
        if (current_line.substr(1, hibf_prefix.size()) == hibf_prefix)
        {
            assert(current_line.substr(hibf_prefix.size() + 2, 11) == "max_bin_id:");
            header.hibf_max_bin_id = current_line.substr(hibf_prefix.size() + 13,
                                                         current_line.size() - hibf_prefix.size() - 13);
        }
    }

    size_t record_idx{};
    do
    {
        auto && record = parse_binning_line(current_line);

        if (starts_with(record.bin_name, merged_bin_prefix))
        {
            // cache low level records in map first to accumulate them
            std::string const prefix(record.bin_name.begin(),
                                     std::find(record.bin_name.begin() + merged_bin_prefix.size() + 1,
                                               record.bin_name.end(),
                                               '_'));
            auto & merged_bin_record = low_level_records[prefix];
            merged_bin_record.filenames.insert(merged_bin_record.filenames.end(),
                                               record.filenames.begin(),
                                               record.filenames.end());

            if (record.max_size > merged_bin_record.max_size)
            {
                merged_bin_record.max_size = record.max_size;
                merged_bin_record.merged_bin_max_size_bin_idx = merged_bin_record.bins;
                merged_bin_record.merged_bin_max_size_filenames = record.filenames;
            }

            merged_bin_record.bins += record.bins;
        }
        else
        {
            // add split record
            records.push_back(record);

            if (record.bin_name == header.hibf_max_bin_id)
                record_idx = records.size() - 1;
        }
    } while (std::getline(binning_file, current_line));

    for (auto & [bin_name, rec] : low_level_records)
    {
        rec.bin_name = bin_name;
        records.push_back(rec);

        if (rec.bin_name.substr(0, header.hibf_max_bin_id.size()) == header.hibf_max_bin_id)
            record_idx = records.size() - 1;
    }

    header.hibf_max_record = &records[record_idx]; // only take a pointer now s.t. it is not invalidated by push_backs
    return std::make_pair(std::move(header), std::move(records));
};
