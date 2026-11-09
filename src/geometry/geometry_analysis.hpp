#pragma once

#include "bounding_box.hpp"
#include "lattice_model.hpp"

#include <cstddef>
#include <optional>
#include <stdexcept>
#include <vector>


namespace ntic::lbm::geometry
{

template <
    lattice::LatticeType Type>
struct GeometryAnalysis2D
{
    using Lattice =
        lattice::Model<Type>;


    static_assert(
        Type ==
            lattice::LatticeType::D2Q9,
        "GeometryAnalysis2D requires D2Q9.");


    double gridSpacing;

    // Optional input domain; otherwise the geometry selects its usual domain.
    std::optional<BoundingBox2D> domainBounds;
    // Analyze the complementary wet region without modifying the geometry.
    bool complemented = false;

    // Actual output domain used for this analysis.
    BoundingBox2D domain;

    std::size_t nx = 0;
    std::size_t ny = 0;

    std::vector<double> scalar;

    std::vector<double> boundaryX;
    std::vector<double> boundaryY;
    
    // Number of geometric boundary edges.
    // Independent of grid resolution.
    std::size_t nBoundaries = 0;
    
    // Boundary entries are grouped by Edge ID.
    // Edge e occupies the half-open interval:
    // [boundaryOffsets[e], boundaryOffsets[e + 1]).
    //
    // Size: nBoundaries + 1.
    // boundaryOffsets[0] == 0.
    // boundaryOffsets.back() == boundaryX.size().
    //
    // A cell may appear in multiple Edge groups.
    // Empty Edge groups are allowed.
    std::vector<std::size_t> boundaryOffsets{0};
    
    std::vector<double> q;


    explicit GeometryAnalysis2D(
        const double gridSpacing)
        :
        gridSpacing(gridSpacing)
    {
        if(gridSpacing <= 0.0)
        {
            throw std::invalid_argument(
                "Grid spacing must be positive.");
        }
    }
};


template <
    lattice::LatticeType Type>
struct GeometryAnalysis3D
{
    using Lattice =
        lattice::Model<Type>;


    static_assert(
        Type ==
            lattice::LatticeType::D3Q15 ||
        Type ==
            lattice::LatticeType::D3Q19 ||
        Type ==
            lattice::LatticeType::D3Q27,
        "GeometryAnalysis3D requires "
        "D3Q15, D3Q19, or D3Q27.");


    double gridSpacing;

    // Optional input domain; otherwise the geometry selects its usual domain.
    std::optional<BoundingBox> domainBounds;
    // Analyze the complementary wet region without modifying the geometry.
    bool complemented = false;

    // Actual output domain used for this analysis.
    BoundingBox domain;

    std::size_t nx = 0;
    std::size_t ny = 0;
    std::size_t nz = 0;

    std::vector<double> scalar;

    std::vector<double> boundaryX;
    std::vector<double> boundaryY;
    std::vector<double> boundaryZ;
    
    // Number of geometric boundary faces.
    // Independent of grid resolution.
    std::size_t nBoundaries = 0;
    
    // Boundary entries are grouped by Face ID.
    // Face f occupies the half-open interval:
    // [boundaryOffsets[f], boundaryOffsets[f + 1]).
    //
    // Size: nBoundaries + 1.
    // boundaryOffsets[0] == 0.
    // boundaryOffsets.back() == boundaryX.size().
    //
    // A cell may appear in multiple Face groups.
    // Empty Face groups are allowed.
    std::vector<std::size_t> boundaryOffsets{0};
    
    std::vector<double> q;


    explicit GeometryAnalysis3D(
        const double gridSpacing)
        :
        gridSpacing(gridSpacing)
    {
        if(gridSpacing <= 0.0)
        {
            throw std::invalid_argument(
                "Grid spacing must be positive.");
        }
    }
};

} // namespace ntic::lbm::geometry
