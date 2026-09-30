#pragma once

#include "bounding_box.hpp"
#include "bounding_box_2d.hpp"
#include "lattice_model.hpp"

#include <cstddef>
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

    BoundingBox2D domain;

    std::size_t nx = 0;
    std::size_t ny = 0;

    std::vector<double> scalar;

    std::vector<double> boundaryX;
    std::vector<double> boundaryY;

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


    std::size_t boundaryCount() const noexcept
    {
        return boundaryX.size();
    }


    void initializeQ()
    {
        if(boundaryX.size() !=
           boundaryY.size())
        {
            throw std::runtime_error(
                "2D boundary coordinate arrays "
                "have inconsistent sizes.");
        }


        q.assign(
            boundaryCount() *
                Lattice::nStencils,
            -1.0);
    }


    double& qValue(
        const std::size_t boundaryID,
        const std::size_t stencilID)
    {
        return q.at(
            boundaryID *
                Lattice::nStencils +
            stencilID);
    }


    const double& qValue(
        const std::size_t boundaryID,
        const std::size_t stencilID) const
    {
        return q.at(
            boundaryID *
                Lattice::nStencils +
            stencilID);
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

    BoundingBox domain;

    std::size_t nx = 0;
    std::size_t ny = 0;
    std::size_t nz = 0;

    std::vector<double> scalar;

    std::vector<double> boundaryX;
    std::vector<double> boundaryY;
    std::vector<double> boundaryZ;

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


    std::size_t boundaryCount() const noexcept
    {
        return boundaryX.size();
    }


    void initializeQ()
    {
        if(boundaryX.size() !=
               boundaryY.size() ||
           boundaryX.size() !=
               boundaryZ.size())
        {
            throw std::runtime_error(
                "3D boundary coordinate arrays "
                "have inconsistent sizes.");
        }


        q.assign(
            boundaryCount() *
                Lattice::nStencils,
            -1.0);
    }


    double& qValue(
        const std::size_t boundaryID,
        const std::size_t stencilID)
    {
        return q.at(
            boundaryID *
                Lattice::nStencils +
            stencilID);
    }


    const double& qValue(
        const std::size_t boundaryID,
        const std::size_t stencilID) const
    {
        return q.at(
            boundaryID *
                Lattice::nStencils +
            stencilID);
    }
};

} // namespace ntic::lbm::geometry