
#pragma once

#include "bounding_box.hpp"
#include "cell_type.hpp"
#include "flow_type.hpp"
#include "geometry_analysis.hpp"
#include "point.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <limits>

namespace ntic::lbm::geometry
{

struct Box
{
    Point min;
    Point max;

    BoundingBox bounds;
    BoundingBox openBox;

    FlowType flowType;


    Box(
        const Point& inputMin,
        const Point& inputMax,
        const FlowType inputFlowType = FlowType::Internal,
        const BoundingBox* inputOpenBox = nullptr)
        : min(inputMin),
          max(inputMax),
          bounds{inputMin, inputMax},
          flowType(inputFlowType)
    {
        if(max[0] <= min[0] ||
           max[1] <= min[1] ||
           max[2] <= min[2])
        {
            throw std::invalid_argument(
                "Box dimensions must be positive.");
        }

        if(flowType == FlowType::External)
        {
            if(inputOpenBox != nullptr)
            {
                validateOpenBox(
                    *inputOpenBox);

                openBox =
                    *inputOpenBox;
            }
            else
            {
                openBox =
                    readOpenBox();
            }
        }
    }


    [[nodiscard]]
    BoundingBox boundingBox() const noexcept
    {
        return bounds;
    }


    void translate(
        const Point* targetPoint = nullptr)
    {
        Point displacement{};

        if(flowType == FlowType::Internal)
        {
            if(targetPoint == nullptr)
            {
                throw std::invalid_argument(
                    "Internal Box translation requires a target point.");
            }

            displacement = {
                (*targetPoint)[0] - bounds.min[0],
                (*targetPoint)[1] - bounds.min[1],
                (*targetPoint)[2] - bounds.min[2]
            };
        }
        else if(flowType == FlowType::External)
        {
            if(targetPoint != nullptr)
            {
                throw std::invalid_argument(
                    "External Box translation does not accept a target point.");
            }

            displacement = {
                -openBox.min[0],
                -openBox.min[1],
                -openBox.min[2]
            };
        }
        else
        {
            throw std::invalid_argument(
                "Unsupported Box flow type.");
        }

        for(std::size_t d = 0;
            d < 3;
            ++d)
        {
            min[d] +=
                displacement[d];

            max[d] +=
                displacement[d];

            bounds.min[d] +=
                displacement[d];

            bounds.max[d] +=
                displacement[d];

            if(flowType ==
               FlowType::External)
            {
                openBox.min[d] +=
                    displacement[d];

                openBox.max[d] +=
                    displacement[d];
            }
        }
    }



    template <lattice::LatticeType Type>
    void analysis(
        GeometryAnalysis3D<Type>& analysis) const
    {
        interiorAreaAnalysis(analysis);
        boundaryQAnalysis(analysis);
    }



    template <lattice::LatticeType Type>
    void boundaryQAnalysis(
        GeometryAnalysis3D<Type>& analysis) const
    {
        using Lattice =
            typename GeometryAnalysis3D<Type>::Lattice;
    
        constexpr std::size_t nLinks =
            Lattice::nStencils - 1;
    
        const std::size_t boundaryCount =
            analysis.boundaryX.size();
    
        if(analysis.boundaryY.size() != boundaryCount ||
           analysis.boundaryZ.size() != boundaryCount)
        {
            throw std::invalid_argument(
                "Boundary coordinate sizes do not match.");
        }

        if(analysis.nBoundaries != 6 || analysis.boundaryOffsets.size() != 7 ||
           analysis.boundaryOffsets.front() != 0 ||
           analysis.boundaryOffsets.back() != boundaryCount ||
           !std::is_sorted(analysis.boundaryOffsets.begin(),
                           analysis.boundaryOffsets.end()))
        {
            throw std::invalid_argument(
                "Box boundary entries must be grouped by Face ID.");
        }
    
        const double h = analysis.gridSpacing;
    
        if(!std::isfinite(h) || h <= 0.0)
        {
            throw std::invalid_argument(
                "Grid spacing must be finite and positive.");
        }
    
        analysis.validateBoundaryCells();

        analysis.q.assign(
            boundaryCount * nLinks, -1.0);
    
        #pragma omp parallel for schedule(static)
        for(std::ptrdiff_t index = 0;
            index < static_cast<std::ptrdiff_t>(boundaryCount);
            ++index)
        {
            const std::size_t b =
                static_cast<std::size_t>(index);
            const auto position=analysis.cellCenter(analysis.boundaryCells[b]);
            const std::size_t faceID = static_cast<std::size_t>(
                std::upper_bound(analysis.boundaryOffsets.begin(),
                                 analysis.boundaryOffsets.end(), b) -
                analysis.boundaryOffsets.begin() - 1);
    
            const double p[3] = {
                position[0],
                position[1],
                position[2]
            };
    
            for(std::size_t i = 1;
                i < Lattice::nStencils;
                ++i)
            {
                const double d[3] = {
                    h * static_cast<double>(Lattice::ex[i]),
                    h * static_cast<double>(Lattice::ey[i]),
                    h * static_cast<double>(Lattice::ez[i])
                };
    
                double tEnter =
                    -std::numeric_limits<double>::infinity();
    
                double tExit =
                    std::numeric_limits<double>::infinity();

                std::size_t enterFace = 6;
                std::size_t exitFace = 6;
    
                bool valid = true;
    
                for(std::size_t axis = 0; axis < 3; ++axis)
                {
                    if(d[axis] == 0.0)
                    {
                        if(p[axis] < min[axis] ||
                           p[axis] > max[axis])
                        {
                            valid = false;
                            break;
                        }
    
                        continue;
                    }
    
                    const double t1 =
                        (min[axis] - p[axis]) / d[axis];
    
                    const double t2 =
                        (max[axis] - p[axis]) / d[axis];
    
                    const double nearT =
                        std::min(t1, t2);
    
                    const double farT =
                        std::max(t1, t2);
    
                    const std::size_t nearFace =
                        2 * axis + (d[axis] > 0.0 ? 0 : 1);
                    const std::size_t farFace =
                        2 * axis + (d[axis] > 0.0 ? 1 : 0);
                    if(nearT > tEnter ||
                       (nearT == tEnter && nearFace < enterFace))
                    {
                        tEnter = nearT;
                        enterFace = nearFace;
                    }
                    if(farT < tExit ||
                       (farT == tExit && farFace < exitFace))
                    {
                        tExit = farT;
                        exitFace = farFace;
                    }
    
                    if(tEnter > tExit)
                    {
                        valid = false;
                        break;
                    }
                }
    
                if(!valid)
                {
                    continue;
                }
    
                double qValue = -1.0;
                std::size_t owningFace = 6;
    
                if(tEnter > 0.0 && tEnter <= 1.0)
                {
                    qValue = tEnter;
                    owningFace = enterFace;
                }
    
                if(tExit > 0.0 && tExit <= 1.0 &&
                   (qValue < 0.0 || tExit < qValue))
                {
                    qValue = tExit;
                    owningFace = exitFace;
                }

                // Exact shared-edge/corner ties use the lowest Face ID.
                if(qValue > 0.0 && tEnter == tExit)
                    owningFace = std::min(enterFace, exitFace);
    
                analysis.q[b * nLinks + (i - 1)] =
                    owningFace == faceID ? qValue : -1.0;
            }
        }
    }



    template <lattice::LatticeType Type>
    void interiorAreaAnalysis(
        GeometryAnalysis3D<Type>& analysis) const
    {
        const bool complemented = analysis.complemented;
        const double h = analysis.gridSpacing;

        if(!std::isfinite(h) || h <= 0.0)
        {
            throw std::invalid_argument(
                "Grid spacing must be finite and positive.");
        }

        analysis.domain = analysis.domainBounds.value_or(
            flowType == FlowType::Internal ? bounds : openBox);

        const BoundingBox& domain = analysis.domain;

        analysis.nx = static_cast<std::size_t>(
            std::ceil(domain.width() / h));

        analysis.ny = static_cast<std::size_t>(
            std::ceil(domain.height() / h));

        analysis.nz = static_cast<std::size_t>(
            std::ceil(domain.depth() / h));

        analysis.validateGridFrame();

        const std::size_t nx = analysis.nx;
        const std::size_t ny = analysis.ny;
        const std::size_t nz = analysis.nz;

        const std::size_t cellCount = nx * ny * nz;

        std::vector<CellType> cellTypes(
            cellCount, CellType::Dry);

        std::vector<std::array<bool, 6>> faceMembership(cellCount);

        // A StencilCell has half-width h, not 0.5*h.
        // Its maximum center-to-corner distance is sqrt(3)*h.
        const double stencilHalfDiagonal =
            std::sqrt(3.0) * h;

        // Classify all cell centers in parallel.
        #pragma omp parallel for schedule(static)
        for(std::ptrdiff_t index = 0;
            index < static_cast<std::ptrdiff_t>(cellCount);
            ++index)
        {
            const std::size_t cellID =
                static_cast<std::size_t>(index);

            const std::size_t i = cellID % nx;
            const std::size_t tmp = cellID / nx;
            const std::size_t j = tmp % ny;
            const std::size_t k = tmp / ny;

            const Point center{
                analysis.centerCoordinate(0,i),
                analysis.centerCoordinate(1,j),
                analysis.centerCoordinate(2,k)
            };

            const bool inside = pointInside(center);

            const bool wet = complemented !=
                (flowType == FlowType::Internal ? inside : !inside);

            // The surface itself is never a wet center,
            // including for External flow.
            if(!wet ||
               distanceToBoundary(center) == 0.0)
            {
                cellTypes[cellID] = CellType::Dry;
                continue;
            }

            // A conservative rejection test:
            // the surface cannot intersect the StencilCell
            // if it is farther than its half-diagonal.
            if(distanceToBoundary(center) >
               stencilHalfDiagonal)
            {
                cellTypes[cellID] = CellType::Interior;
                continue;
            }

            // Closed StencilCell: center +/- h.
            const Point stencilMin{
                center[0] - h,
                center[1] - h,
                center[2] - h
            };

            const Point stencilMax{
                center[0] + h,
                center[1] + h,
                center[2] + h
            };

            faceMembership[cellID] = boundaryIntersections(stencilMin, stencilMax);
            const auto& membership = faceMembership[cellID];
            const bool intersects = std::any_of(
                membership.begin(), membership.end(), [](bool hit) { return hit; });
            
            cellTypes[cellID] =
                intersects ? CellType::Boundary : CellType::Interior;
        }

        // Face IDs: XMin, XMax, YMin, YMax, ZMin, ZMax.
        analysis.nBoundaries = 6;
        
        analysis.boundaryOffsets.assign(
            analysis.nBoundaries + 1, 0);
        
        // Histogram.
        for(const auto& membership : faceMembership)
        {
            for(std::size_t face = 0; face < analysis.nBoundaries; ++face)
                if(membership[face]) ++analysis.boundaryOffsets[face + 1];
        }
        
        // Prefix sum.
        for(std::size_t k = 0; k < analysis.nBoundaries; ++k)
        {
            analysis.boundaryOffsets[k + 1] +=
                analysis.boundaryOffsets[k];
        }
        
        const std::size_t boundaryCount = analysis.boundaryOffsets.back();
        analysis.boundaryCells.resize(boundaryCount);
        analysis.boundaryX.resize(boundaryCount);
        analysis.boundaryY.resize(boundaryCount);
        analysis.boundaryZ.resize(boundaryCount);

        // Disjoint Face ranges, each in ascending physical cell-index order.
        #pragma omp parallel for schedule(static)
        for(int face = 0; face < 6; ++face)
        {
            std::size_t slot = analysis.boundaryOffsets[face];
            for(std::size_t id = 0; id < cellCount; ++id)
            {
                if(!faceMembership[id][face]) continue;
                const std::size_t i = id % nx;
                const std::size_t j = (id / nx) % ny;
                const std::size_t k = id / (nx * ny);
                analysis.boundaryCells[slot]=id;
                analysis.boundaryX[slot] = analysis.centerCoordinate(0,i);
                analysis.boundaryY[slot] = analysis.centerCoordinate(1,j);
                analysis.boundaryZ[slot] = analysis.centerCoordinate(2,k);
                ++slot;
            }
        }

        // Keep the existing scalar representation
        // for compatibility with VTK output.
        analysis.scalar.resize(cellCount);

        #pragma omp parallel for schedule(static)
        for(std::ptrdiff_t index = 0;
            index < static_cast<std::ptrdiff_t>(cellCount);
            ++index)
        {
            const std::size_t cellID =
                static_cast<std::size_t>(index);

            analysis.scalar[cellID] =
                static_cast<double>(
                    static_cast<int>(cellTypes[cellID]));
        }
    }


private:


    [[nodiscard]]
    bool pointInside(
        const Point& point) const noexcept
    {
        return
            point[0] > min[0] &&
            point[0] < max[0] &&
            point[1] > min[1] &&
            point[1] < max[1] &&
            point[2] > min[2] &&
            point[2] < max[2];
    }


    [[nodiscard]]
    double distanceToBoundary(
        const Point& point) const noexcept
    {
        const double dx =
            std::max(
                std::max(
                    min[0] - point[0],
                    0.0),
                point[0] - max[0]);

        const double dy =
            std::max(
                std::max(
                    min[1] - point[1],
                    0.0),
                point[1] - max[1]);

        const double dz =
            std::max(
                std::max(
                    min[2] - point[2],
                    0.0),
                point[2] - max[2]);

        if(dx > 0.0 ||
           dy > 0.0 ||
           dz > 0.0)
        {
            return
                std::sqrt(
                    dx * dx +
                    dy * dy +
                    dz * dz);
        }

        return
            std::min({
                point[0] - min[0],
                max[0] - point[0],
                point[1] - min[1],
                max[1] - point[1],
                point[2] - min[2],
                max[2] - point[2]
            });
    }


    // Stable Face IDs: 0 XMin, 1 XMax, 2 YMin, 3 YMax, 4 ZMin, 5 ZMax.
    // Closed contact at edges and vertices belongs to all incident faces.
    [[nodiscard]]
    std::array<bool, 6> boundaryIntersections(
        const Point& cellMin, const Point& cellMax) const noexcept
    {
        std::array<bool, 6> result{};
        for(std::size_t face = 0; face < result.size(); ++face)
        {
            const std::size_t axis = face / 2;
            const double plane = face % 2 == 0 ? min[axis] : max[axis];
            bool hit = cellMin[axis] <= plane && cellMax[axis] >= plane;
            for(std::size_t d = 0; d < 3; ++d)
            {
                if(d != axis && (cellMax[d] < min[d] || cellMin[d] > max[d]))
                    hit = false;
            }
            result[face] = hit;
        }
        return result;
    }

    void validateOpenBox(
        const BoundingBox& box) const
    {
        if(box.width() <= 0.0 ||
           box.height() <= 0.0 ||
           box.depth() <= 0.0)
        {
            throw std::invalid_argument(
                "Open box dimensions must be positive.");
        }

        if(!box.contains(bounds))
        {
            throw std::invalid_argument(
                "Open box must contain the Box bounds.");
        }
    }


    [[nodiscard]]
    BoundingBox readOpenBox() const
    {
        std::cout
            << "Box bounds:\n"
            << "  min = ("
            << bounds.min[0] << ", "
            << bounds.min[1] << ", "
            << bounds.min[2] << ")\n"
            << "  max = ("
            << bounds.max[0] << ", "
            << bounds.max[1] << ", "
            << bounds.max[2] << ")\n";

        Point origin{};
        Point size{};

        std::cout
            << "Enter open box origin x y z: ";

        if(!(std::cin >>
             origin[0] >>
             origin[1] >>
             origin[2]))
        {
            throw std::runtime_error(
                "Failed to read open box origin.");
        }

        std::cout
            << "Enter open box size x y z: ";

        if(!(std::cin >>
             size[0] >>
             size[1] >>
             size[2]))
        {
            throw std::runtime_error(
                "Failed to read open box size.");
        }

        const BoundingBox box{
            origin,
            {
                origin[0] + size[0],
                origin[1] + size[1],
                origin[2] + size[2]
            }
        };

        validateOpenBox(box);

        return box;
    }
};

} // namespace ntic::lbm::geometry
