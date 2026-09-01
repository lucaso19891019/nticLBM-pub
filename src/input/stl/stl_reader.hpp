#pragma once

#include <array>
#include <cstddef>
#include <filesystem>
#include <vector>

namespace ntic::lbm::stl {

//=============================================================================
// STL types
//=============================================================================

enum class STLFormat
{
    ASCII,
    Binary
};

using STLVector = std::array<double, 3>;

struct STLFacet
{
    STLVector normal{};
    std::array<STLVector, 3> vertices{};
};

struct STLData
{
    STLFormat format = STLFormat::ASCII;
    std::vector<STLFacet> facets;

    [[nodiscard]]
    std::size_t size() const noexcept
    {
        return facets.size();
    }

    [[nodiscard]]
    bool empty() const noexcept
    {
        return facets.empty();
    }
};


//=============================================================================
// STL reader
//=============================================================================

[[nodiscard]]
STLData read(const std::filesystem::path& path);

} // namespace ntic::lbm::stl