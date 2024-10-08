#pragma once

#include <filesystem>
#include <iostream>

inline void print_debug_file(std::filesystem::path const & path, size_t const column_width = 18u)
{
    std::ios::fmtflags const original_flags = std::cout.setf(std::ios::left, std::ios::adjustfield);

    std::ifstream layouting_file{path};

    std::string buffer{};
    buffer.reserve(column_width);
    auto flush_buffer = [&buffer, &column_width]()
    {
        std::cout.width(column_width);
        std::cout << buffer;
        buffer.clear();
    };

    for (auto it = std::istreambuf_iterator<char>(layouting_file); it != std::istreambuf_iterator<char>(); ++it)
    {
        char const letter{*it};
        switch (letter)
        {
        case '\n':
            flush_buffer();
            std::cout << '\n';
            break;
        case '\t':
            flush_buffer();
            break;
        default:
            buffer += letter;
            break;
        }
    }

    std::cout.setf(original_flags);
}
