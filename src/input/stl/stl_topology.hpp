#pragma once

#include <cstddef>
#include <vector>

namespace ntic::lbm::stl
{

struct GeometryBounds
{
    std::array<double,3> min;
    std::array<double,3> max;
    double scale;
};


struct STLComponent
{
    std::vector<std::size_t> facets;

    GeometryBounds bounds;

    double signedVolume = 0.0;
};


using STLComponents =
    std::vector<STLComponent>;


struct EdgeUse
{
    std::size_t facetID;

    bool forward;
};


struct EdgeKey
{
    std::size_t v0;
    std::size_t v1;

    bool operator==(const EdgeKey& other) const noexcept
    {
        return v0 == other.v0 &&
               v1 == other.v1;
    }
};


struct FacetTopology
{
    std::vector<std::vector<std::size_t>>
        adjacency;

    STLComponents components;
};

}