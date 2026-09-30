#pragma once

#include "bounding_box.hpp"

#include <cstddef>
#include <stdexcept>
#include <vector>


namespace ntic::lbm::geometry
{

//=============================================================================
// GeometryAnalysis2D
//=============================================================================

template <typename LatticeModel>
struct GeometryAnalysis2D
{
    double gridSpacing;

    BoundingBox2D domain;

    std::size_t nx =
        0;

    std::size_t ny =
        0;

    std::vector<double> scalar;

    std::vector<double> boundaryX;

    std::vector<double> boundaryY;

    // Lattice-model-dependent boundary-link q data will be added here.


    explicit GeometryAnalysis2D(
        const double gridSpacing)
        : gridSpacing(gridSpacing)
    {
        if(gridSpacing <= 0.0)
        {
            throw std::invalid_argument(
                "Grid spacing must be positive.");
        }
    }
};


//=============================================================================
// GeometryAnalysis3D
//=============================================================================

template <typename LatticeModel>
struct GeometryAnalysis3D
{
    double gridSpacing;

    BoundingBox domain;

    std::size_t nx =
        0;

    std::size_t ny =
        0;

    std::size_t nz =
        0;

    std::vector<double> scalar;

    std::vector<double> boundaryX;

    std::vector<double> boundaryY;

    std::vector<double> boundaryZ;

    // Lattice-model-dependent boundary-link q data will be added here.


    explicit GeometryAnalysis3D(
        const double gridSpacing)
        : gridSpacing(gridSpacing)
    {
        if(gridSpacing <= 0.0)
        {
            throw std::invalid_argument(
                "Grid spacing must be positive.");
        }
    }
};

} // namespace ntic::lbm::geometry
