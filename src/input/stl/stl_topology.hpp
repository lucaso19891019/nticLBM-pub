#pragma once

#include <array>
#include <cstddef>
#include <vector>

namespace ntic::lbm::stl
{

struct GeometricVertices
{
    std::vector<std::array<double,3>>
        vertices;

    std::vector<std::array<std::size_t,3>>
        facetVertexIDs;
};


struct GeometryBounds
{
    std::array<double,3> min;
    std::array<double,3> max;
    std::array<double,3> center;

    double scale;
};


struct STLComponent
{
    std::vector<std::size_t>
        facets;

    GeometryBounds bounds;

    double signedVolume = 0.0;
};


struct STLEdge
{
    std::size_t v0;

    std::size_t v1;

    std::array<std::size_t,2>
        facets;
};


struct FacetGeometry
{
    std::array<double,3> centroid{};

    std::array<double,3> normal{};

    double area = 0.0;

    double quality = 0.0;

    double centroidRadius = 0.0;
};


struct FacetTopology
{
    GeometricVertices geometry;


    std::vector<FacetGeometry>
        facetGeometry;


    std::vector<std::array<std::size_t,3>>
        adjacency;


    std::vector<STLEdge>
        edges;


    std::vector<std::size_t>
        facetComponentIDs;


    STLComponents components;
};

}