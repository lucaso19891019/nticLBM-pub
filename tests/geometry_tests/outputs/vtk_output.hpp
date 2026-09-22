#pragma once

#include <filesystem>
#include <stdexcept>
#include <system_error>


namespace ntic::lbm::geometry::test
{

inline void recreateOutputDirectory(
    const std::filesystem::path& path)
{
    std::error_code error;


    const bool exists =
        std::filesystem::exists(
            path,
            error);

    if(error)
    {
        throw std::runtime_error(
            "Failed to check output directory: " +
            path.string());
    }


    if(exists)
    {
        std::filesystem::remove_all(
            path,
            error);

        if(error)
        {
            throw std::runtime_error(
                "Failed to remove output directory: " +
                path.string());
        }
    }


    std::filesystem::create_directories(
        path,
        error);

    if(error)
    {
        throw std::runtime_error(
            "Failed to create output directory: " +
            path.string());
    }
}

} // namespace ntic::lbm::geometry::test