
#pragma once

#include "bounding_box.hpp"
#include "cell_type.hpp"
#include "flow_type.hpp"
#include "geometry_analysis.hpp"
#include "point.hpp"

#include <algorithm>
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
    
        const double h = analysis.gridSpacing;
    
        if(!std::isfinite(h) || h <= 0.0)
        {
            throw std::invalid_argument(
                "Grid spacing must be finite and positive.");
        }
    
        analysis.q.assign(
            boundaryCount * nLinks, -1.0);
    
        #pragma omp parallel for schedule(static)
        for(std::ptrdiff_t index = 0;
            index < static_cast<std::ptrdiff_t>(boundaryCount);
            ++index)
        {
            const std::size_t b =
                static_cast<std::size_t>(index);
    
            const double p[3] = {
                analysis.boundaryX[b],
                analysis.boundaryY[b],
                analysis.boundaryZ[b]
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
    
                    tEnter =
                        std::max(tEnter, nearT);
    
                    tExit =
                        std::min(tExit, farT);
    
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
    
                if(tEnter > 0.0 && tEnter <= 1.0)
                {
                    qValue = tEnter;
                }
    
                if(tExit > 0.0 && tExit <= 1.0 &&
                   (qValue < 0.0 || tExit < qValue))
                {
                    qValue = tExit;
                }
    
                analysis.q[b * nLinks + (i - 1)] =
                    qValue;
            }
        }
    }



    template <lattice::LatticeType Type>
    void interiorAreaAnalysis(
        GeometryAnalysis3D<Type>& analysis) const
    {
        const double h = analysis.gridSpacing;

        if(!std::isfinite(h) || h <= 0.0)
        {
            throw std::invalid_argument(
                "Grid spacing must be finite and positive.");
        }

        analysis.domain =
            flowType == FlowType::Internal
                ? bounds
                : openBox;

        const BoundingBox& domain = analysis.domain;

        analysis.nx = static_cast<std::size_t>(
            std::ceil(domain.width() / h));

        analysis.ny = static_cast<std::size_t>(
            std::ceil(domain.height() / h));

        analysis.nz = static_cast<std::size_t>(
            std::ceil(domain.depth() / h));

        const std::size_t nx = analysis.nx;
        const std::size_t ny = analysis.ny;
        const std::size_t nz = analysis.nz;

        const std::size_t cellCount = nx * ny * nz;

        std::vector<CellType> cellTypes(
            cellCount, CellType::Dry);

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
                domain.min[0] +
                    (static_cast<double>(i) + 0.5) * h,
                domain.min[1] +
                    (static_cast<double>(j) + 0.5) * h,
                domain.min[2] +
                    (static_cast<double>(k) + 0.5) * h
            };

            const bool inside = pointInside(center);

            const bool wet =
                flowType == FlowType::Internal
                    ? inside
                    : !inside;

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

            const SurfaceRelation relation =
                surfaceRelation(stencilMin, stencilMax);

            cellTypes[cellID] =
                relation == SurfaceRelation::None
                    ? CellType::Interior
                    : CellType::Boundary;
        }

        // Compact Boundary coordinates in parallel.
        // Fixed-size blocks avoid atomic contention and
        // preserve a deterministic cell-index order.
        constexpr std::size_t blockSize = 4096;

        const std::size_t blockCount =
            cellCount / blockSize +
            (cellCount % blockSize != 0 ? 1 : 0);

        std::vector<std::size_t> blockOffsets(
            blockCount + 1, 0);

        #pragma omp parallel for schedule(static)
        for(std::ptrdiff_t block = 0;
            block < static_cast<std::ptrdiff_t>(blockCount);
            ++block)
        {
            const std::size_t b =
                static_cast<std::size_t>(block);

            const std::size_t begin = b * blockSize;
            const std::size_t end =
                std::min(begin + blockSize, cellCount);

            std::size_t count = 0;

            for(std::size_t cellID = begin;
                cellID < end;
                ++cellID)
            {
                if(cellTypes[cellID] == CellType::Boundary)
                {
                    ++count;
                }
            }

            blockOffsets[b + 1] = count;
        }

        // Prefix sum: one small serial pass over blocks.
        for(std::size_t b = 0;
            b < blockCount;
            ++b)
        {
            blockOffsets[b + 1] += blockOffsets[b];
        }

        const std::size_t boundaryCount =
            blockOffsets[blockCount];

        analysis.boundaryX.resize(boundaryCount);
        analysis.boundaryY.resize(boundaryCount);
        analysis.boundaryZ.resize(boundaryCount);

        #pragma omp parallel for schedule(static)
        for(std::ptrdiff_t block = 0;
            block < static_cast<std::ptrdiff_t>(blockCount);
            ++block)
        {
            const std::size_t b =
                static_cast<std::size_t>(block);

            const std::size_t begin = b * blockSize;
            const std::size_t end =
                std::min(begin + blockSize, cellCount);

            std::size_t slot = blockOffsets[b];

            for(std::size_t cellID = begin;
                cellID < end;
                ++cellID)
            {
                if(cellTypes[cellID] != CellType::Boundary)
                {
                    continue;
                }

                const std::size_t i = cellID % nx;
                const std::size_t tmp = cellID / nx;
                const std::size_t j = tmp % ny;
                const std::size_t k = tmp / ny;

                analysis.boundaryX[slot] =
                    domain.min[0] +
                    (static_cast<double>(i) + 0.5) * h;

                analysis.boundaryY[slot] =
                    domain.min[1] +
                    (static_cast<double>(j) + 0.5) * h;

                analysis.boundaryZ[slot] =
                    domain.min[2] +
                    (static_cast<double>(k) + 0.5) * h;

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

    enum class SurfaceRelation
    {
        None,
        Touch,
        Cross
    };


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


    [[nodiscard]]
    SurfaceRelation surfaceRelation(
        const Point& cellMin,
        const Point& cellMax) const noexcept
    {
        const bool closedOverlap =
            cellMax[0] >= min[0] &&
            cellMin[0] <= max[0] &&
            cellMax[1] >= min[1] &&
            cellMin[1] <= max[1] &&
            cellMax[2] >= min[2] &&
            cellMin[2] <= max[2];

        if(!closedOverlap)
        {
            return
                SurfaceRelation::None;
        }

        const bool cellInsideOpen =
            cellMin[0] > min[0] &&
            cellMax[0] < max[0] &&
            cellMin[1] > min[1] &&
            cellMax[1] < max[1] &&
            cellMin[2] > min[2] &&
            cellMax[2] < max[2];

        if(cellInsideOpen)
        {
            return
                SurfaceRelation::None;
        }

        const bool cellInsideClosed =
            cellMin[0] >= min[0] &&
            cellMax[0] <= max[0] &&
            cellMin[1] >= min[1] &&
            cellMax[1] <= max[1] &&
            cellMin[2] >= min[2] &&
            cellMax[2] <= max[2];

        const bool openOverlap =
            cellMax[0] > min[0] &&
            cellMin[0] < max[0] &&
            cellMax[1] > min[1] &&
            cellMin[1] < max[1] &&
            cellMax[2] > min[2] &&
            cellMin[2] < max[2];

        if(openOverlap &&
           !cellInsideClosed)
        {
            return
                SurfaceRelation::Cross;
        }

        return
            SurfaceRelation::Touch;
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
