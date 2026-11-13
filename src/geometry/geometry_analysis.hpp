#pragma once

#include "bounding_box.hpp"
#include "cell_type.hpp"
#include "lattice_model.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <stdexcept>
#include <vector>


namespace ntic::lbm::geometry
{

namespace detail
{
// One explicit arithmetic path for global and tile-local cell centers.
inline double gridCenterCoordinate(double origin,std::size_t index,double spacing) noexcept
{ return std::fma(static_cast<double>(index)+0.5,spacing,origin); }
} // namespace detail

template<std::size_t D>
struct AnalysisGridFrame
{
    std::array<double,D> origin{};
    std::array<std::size_t,D> offset{};

    void validate(const std::array<std::size_t,D>& dimensions,double spacing) const
    {
        if(!(spacing>0) || !std::isfinite(spacing))
            throw std::invalid_argument("Invalid anchored grid spacing.");
        for(std::size_t d=0; d<D; ++d)
        {
            if(!std::isfinite(origin[d]) || dimensions[d]==0 ||
               offset[d]>std::numeric_limits<std::size_t>::max()-(dimensions[d]-1))
                throw std::invalid_argument("Invalid anchored grid origin or index range.");
            const auto last=offset[d]+dimensions[d]-1;
            // Preserve an exactly representable half-integer before fma.
            if(static_cast<double>(last)>=std::ldexp(1.0,std::numeric_limits<double>::digits-1) ||
               !std::isfinite(detail::gridCenterCoordinate(origin[d],offset[d],spacing)) ||
               !std::isfinite(detail::gridCenterCoordinate(origin[d],last,spacing)))
                throw std::invalid_argument("Anchored grid coordinates are not representable.");
        }
    }
};

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
    // Optional global lattice anchor; offsets address the first cell of this window.
    std::optional<AnalysisGridFrame<2>> gridFrame;
    // Analyze the complementary wet region without modifying the geometry.
    bool complemented = false;

    // Actual output domain used for this analysis.
    BoundingBox2D domain;

    std::size_t nx = 0;
    std::size_t ny = 0;

    std::vector<double> scalar;

    // Window-local scalar indices, one per Face/edge record; duplicates
    // across different groups refer to the same physical cell.
    std::vector<std::size_t> boundaryCells;
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
    // boundaryOffsets.back() == boundaryCells.size() == boundaryX.size().
    //
    // A cell may appear in multiple Edge groups.
    // Empty Edge groups are allowed.
    std::vector<std::size_t> boundaryOffsets{0};
    
    std::vector<double> q;


    double centerCoordinate(std::size_t axis,std::size_t index) const noexcept
    {
        return gridFrame ? detail::gridCenterCoordinate(gridFrame->origin[axis],gridFrame->offset[axis]+index,gridSpacing) :
                           domain.min[axis]+(static_cast<double>(index)+0.5)*gridSpacing;
    }

    void validateGridFrame() const
    {
        if(gridFrame) gridFrame->validate({nx,ny},gridSpacing);
    }


    // Requires validated nonzero dimensions and an in-range scalar index.
    std::array<double,2> cellCenter(std::size_t id) const noexcept
    { return {centerCoordinate(0,id%nx),centerCoordinate(1,(id/nx)%ny)}; }

    void validateBoundaryCells() const
    {
        if(!(gridSpacing>0) || !std::isfinite(gridSpacing) || nx==0 || ny==0 ||
           nx>std::numeric_limits<std::size_t>::max()/ny || scalar.size()!=nx*ny)
            throw std::invalid_argument("Boundary indices require a valid scalar grid.");
        validateGridFrame();
        if(boundaryCells.size()!=boundaryX.size() || boundaryY.size()!=boundaryCells.size())
            throw std::invalid_argument("Boundary index and coordinate counts do not match.");
        for(std::size_t entry=0; entry<boundaryCells.size(); ++entry)
        {
            const auto id=boundaryCells[entry];
            if(id>=scalar.size() || scalar[id]!=static_cast<double>(CellType::Boundary))
                throw std::invalid_argument("Boundary index does not identify a Boundary cell.");
            const auto point=cellCenter(id);
            for(double coordinate:point)
                if(!std::isfinite(coordinate)) throw std::invalid_argument("Non-finite boundary center.");
            if(boundaryX[entry]!=point[0] || boundaryY[entry]!=point[1])
                throw std::invalid_argument("Boundary coordinate does not match its cell index.");
        }
    }


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
    // Optional global lattice anchor; offsets address the first cell of this window.
    std::optional<AnalysisGridFrame<3>> gridFrame;
    // Analyze the complementary wet region without modifying the geometry.
    bool complemented = false;

    // Actual output domain used for this analysis.
    BoundingBox domain;

    std::size_t nx = 0;
    std::size_t ny = 0;
    std::size_t nz = 0;

    std::vector<double> scalar;

    // Window-local scalar indices, one per Face/edge record; duplicates
    // across different groups refer to the same physical cell.
    std::vector<std::size_t> boundaryCells;
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
    // boundaryOffsets.back() == boundaryCells.size() == boundaryX.size().
    //
    // A cell may appear in multiple Face groups.
    // Empty Face groups are allowed.
    std::vector<std::size_t> boundaryOffsets{0};
    
    std::vector<double> q;


    double centerCoordinate(std::size_t axis,std::size_t index) const noexcept
    {
        return gridFrame ? detail::gridCenterCoordinate(gridFrame->origin[axis],gridFrame->offset[axis]+index,gridSpacing) :
                           domain.min[axis]+(static_cast<double>(index)+0.5)*gridSpacing;
    }

    void validateGridFrame() const
    {
        if(gridFrame) gridFrame->validate({nx,ny,nz},gridSpacing);
    }


    // Requires validated nonzero dimensions and an in-range scalar index.
    std::array<double,3> cellCenter(std::size_t id) const noexcept
    { return {centerCoordinate(0,id%nx),centerCoordinate(1,(id/nx)%ny),centerCoordinate(2,id/(nx*ny))}; }

    void validateBoundaryCells() const
    {
        if(!(gridSpacing>0) || !std::isfinite(gridSpacing) || nx==0 || ny==0 || nz==0 ||
           nx>std::numeric_limits<std::size_t>::max()/ny || nx*ny>std::numeric_limits<std::size_t>::max()/nz || scalar.size()!=nx*ny*nz)
            throw std::invalid_argument("Boundary indices require a valid scalar grid.");
        validateGridFrame();
        if(boundaryCells.size()!=boundaryX.size() || boundaryY.size()!=boundaryCells.size() || boundaryZ.size()!=boundaryCells.size())
            throw std::invalid_argument("Boundary index and coordinate counts do not match.");
        for(std::size_t entry=0; entry<boundaryCells.size(); ++entry)
        {
            const auto id=boundaryCells[entry];
            if(id>=scalar.size() || scalar[id]!=static_cast<double>(CellType::Boundary))
                throw std::invalid_argument("Boundary index does not identify a Boundary cell.");
            const auto point=cellCenter(id);
            for(double coordinate:point)
                if(!std::isfinite(coordinate)) throw std::invalid_argument("Non-finite boundary center.");
            if(boundaryX[entry]!=point[0] || boundaryY[entry]!=point[1] || boundaryZ[entry]!=point[2])
                throw std::invalid_argument("Boundary coordinate does not match its cell index.");
        }
    }


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
