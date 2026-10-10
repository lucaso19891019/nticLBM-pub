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

namespace ntic::lbm::geometry
{

struct Sphere
{
    Point center;
    double radius;

    BoundingBox bounds;
    BoundingBox openBox;

    FlowType flowType;


    Sphere(
        const Point& inputCenter,
        const double inputRadius,
        const FlowType inputFlowType = FlowType::Internal,
        const BoundingBox* inputOpenBox = nullptr)
        : center(inputCenter),
          radius(inputRadius),
          bounds{
              {
                  inputCenter[0] - inputRadius,
                  inputCenter[1] - inputRadius,
                  inputCenter[2] - inputRadius
              },
              {
                  inputCenter[0] + inputRadius,
                  inputCenter[1] + inputRadius,
                  inputCenter[2] + inputRadius
              }
          },
          flowType(inputFlowType)
    {
        if(radius <= 0.0)
        {
            throw std::invalid_argument(
                "Sphere radius must be positive.");
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
                    "Internal Sphere translation requires a target point.");
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
                    "External Sphere translation does not accept a target point.");
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
                "Unsupported Sphere flow type.");
        }

        for(std::size_t d = 0;
            d < 3;
            ++d)
        {
            center[d] +=
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
    void analysis(GeometryAnalysis3D<Type>& analysis) const
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
    
            const double px =
                analysis.boundaryX[b] - center[0];
    
            const double py =
                analysis.boundaryY[b] - center[1];
    
            const double pz =
                analysis.boundaryZ[b] - center[2];
    
            const double C =
                px * px + py * py + pz * pz -
                radius * radius;
    
            for(std::size_t i = 1;
                i < Lattice::nStencils;
                ++i)
            {
                const double dx =
                    h * static_cast<double>(Lattice::ex[i]);
    
                const double dy =
                    h * static_cast<double>(Lattice::ey[i]);
    
                const double dz =
                    h * static_cast<double>(Lattice::ez[i]);
    
                const double A =
                    dx * dx + dy * dy + dz * dz;
    
                const double B =
                    2.0 * (px * dx + py * dy + pz * dz);
    
                const double discriminant =
                    B * B - 4.0 * A * C;
    
                if(discriminant < 0.0)
                {
                    continue;
                }
    
                const double sqrtD =
                    std::sqrt(discriminant);
    
                const double rootTerm =
                    -0.5 * (B + std::copysign(sqrtD, B));
    
                double qValue = -1.0;
    
                if(rootTerm == 0.0)
                {
                    const double root =
                        -B / (2.0 * A);
    
                    if(root > 0.0 && root <= 1.0)
                    {
                        qValue = root;
                    }
                }
                else
                {
                    const double root1 =
                        rootTerm / A;
    
                    const double root2 =
                        C / rootTerm;
    
                    if(root1 > 0.0 && root1 <= 1.0)
                    {
                        qValue = root1;
                    }
    
                    if(root2 > 0.0 && root2 <= 1.0 &&
                       (qValue < 0.0 || root2 < qValue))
                    {
                        qValue = root2;
                    }
                }
    
                analysis.q[b * nLinks + (i - 1)] =
                    qValue;
            }
        }
    }


    template <lattice::LatticeType Type>
    void interiorAreaAnalysis(GeometryAnalysis3D<Type>& analysis) const
    {
        const double h = analysis.gridSpacing;
        if(!std::isfinite(h) || h <= 0.0)
        {
            throw std::invalid_argument("Grid spacing must be finite and positive.");
        }

        analysis.domain = flowType == FlowType::Internal ? bounds : openBox;
        const BoundingBox& domain = analysis.domain;
        analysis.nx = static_cast<std::size_t>(std::ceil(domain.width() / h));
        analysis.ny = static_cast<std::size_t>(std::ceil(domain.height() / h));
        analysis.nz = static_cast<std::size_t>(std::ceil(domain.depth() / h));

        const std::size_t nx = analysis.nx;
        const std::size_t ny = analysis.ny;
        const std::size_t nz = analysis.nz;
        const std::size_t cellCount = nx * ny * nz;
        std::vector<CellType> cellTypes(cellCount, CellType::Dry);

        // StencilCell half-width is h; its half-diagonal is sqrt(3)*h.
        const double stencilHalfDiagonal = std::sqrt(3.0) * h;

        #pragma omp parallel for schedule(static)
        for(std::ptrdiff_t index = 0;
            index < static_cast<std::ptrdiff_t>(cellCount); ++index)
        {
            const std::size_t id = static_cast<std::size_t>(index);
            const std::size_t i = id % nx;
            const std::size_t j = (id / nx) % ny;
            const std::size_t k = id / (nx * ny);
            const Point centerPoint{
                domain.min[0] + (static_cast<double>(i) + 0.5) * h,
                domain.min[1] + (static_cast<double>(j) + 0.5) * h,
                domain.min[2] + (static_cast<double>(k) + 0.5) * h
            };

            const bool inside = pointInside(centerPoint);
            const bool wet = flowType == FlowType::Internal ? inside : !inside;
            const double surfaceDistance = distanceToBoundary(centerPoint);

            // The center must be strictly wet and must not lie on the surface.
            if(!wet || surfaceDistance == 0.0)
            {
                cellTypes[id] = CellType::Dry;
                continue;
            }
            if(surfaceDistance > stencilHalfDiagonal)
            {
                cellTypes[id] = CellType::Interior;
                continue;
            }

            const Point stencilMin{
                centerPoint[0] - h, centerPoint[1] - h, centerPoint[2] - h
            };
            const Point stencilMax{
                centerPoint[0] + h, centerPoint[1] + h, centerPoint[2] + h
            };
            const SurfaceRelation relation = surfaceRelation(stencilMin, stencilMax);
            cellTypes[id] = relation == SurfaceRelation::None
                ? CellType::Interior : CellType::Boundary;
        }

        // Deterministic, parallel boundary-coordinate compaction.
        constexpr std::size_t blockSize = 4096;
        const std::size_t blockCount =
            cellCount / blockSize + (cellCount % blockSize != 0 ? 1 : 0);
        std::vector<std::size_t> offsets(blockCount + 1, 0);

        #pragma omp parallel for schedule(static)
        for(std::ptrdiff_t block = 0;
            block < static_cast<std::ptrdiff_t>(blockCount); ++block)
        {
            const std::size_t b = static_cast<std::size_t>(block);
            const std::size_t begin = b * blockSize;
            const std::size_t end = std::min(begin + blockSize, cellCount);
            std::size_t count = 0;
            for(std::size_t id = begin; id < end; ++id)
            {
                if(cellTypes[id] == CellType::Boundary) ++count;
            }
            offsets[b + 1] = count;
        }
        for(std::size_t b = 0; b < blockCount; ++b)
        {
            offsets[b + 1] += offsets[b];
        }

        const std::size_t boundaryCount = offsets[blockCount];
        analysis.boundaryX.resize(boundaryCount);
        analysis.boundaryY.resize(boundaryCount);
        analysis.boundaryZ.resize(boundaryCount);

        #pragma omp parallel for schedule(static)
        for(std::ptrdiff_t block = 0;
            block < static_cast<std::ptrdiff_t>(blockCount); ++block)
        {
            const std::size_t b = static_cast<std::size_t>(block);
            const std::size_t begin = b * blockSize;
            const std::size_t end = std::min(begin + blockSize, cellCount);
            std::size_t slot = offsets[b];
            for(std::size_t id = begin; id < end; ++id)
            {
                if(cellTypes[id] != CellType::Boundary) continue;
                const std::size_t i = id % nx;
                const std::size_t j = (id / nx) % ny;
                const std::size_t k = id / (nx * ny);
                analysis.boundaryX[slot] =
                    domain.min[0] + (static_cast<double>(i) + 0.5) * h;
                analysis.boundaryY[slot] =
                    domain.min[1] + (static_cast<double>(j) + 0.5) * h;
                analysis.boundaryZ[slot] =
                    domain.min[2] + (static_cast<double>(k) + 0.5) * h;
                ++slot;
            }
        }

        // Boundary element indices:
        // 0: Spherical surface.
        //
        // boundaryOffsets[0]: Start of spherical surface.
        // boundaryOffsets[1]: End of spherical surface.
        analysis.nBoundaries = 1;
        analysis.boundaryOffsets = {
            0,
            boundaryCount
        };

        analysis.scalar.resize(cellCount);
        #pragma omp parallel for schedule(static)
        for(std::ptrdiff_t index = 0;
            index < static_cast<std::ptrdiff_t>(cellCount); ++index)
        {
            const std::size_t id = static_cast<std::size_t>(index);
            analysis.scalar[id] = static_cast<double>(static_cast<int>(cellTypes[id]));
        }
    }


private:


    [[nodiscard]]
    bool pointInside(
        const Point& point) const noexcept
    {
        const double dx =
            point[0] - center[0];

        const double dy =
            point[1] - center[1];

        const double dz =
            point[2] - center[2];

        return
            dx * dx +
            dy * dy +
            dz * dz <
            radius * radius;
    }


    [[nodiscard]]
    double distanceToBoundary(
        const Point& point) const noexcept
    {
        const double dx =
            point[0] - center[0];

        const double dy =
            point[1] - center[1];

        const double dz =
            point[2] - center[2];

        const double distance =
            std::sqrt(
                dx * dx +
                dy * dy +
                dz * dz);

        return
            std::abs(
                distance - radius);
    }


    [[nodiscard]]
    SurfaceRelation surfaceRelation(
        const Point& cellMin,
        const Point& cellMax) const noexcept
    {
        const double radiusSquared = radius * radius;
        double minDistanceSquared =
            0.0;

        double maxDistanceSquared =
            0.0;

        for(std::size_t d = 0;
            d < 3;
            ++d)
        {
            const double closest =
                std::clamp(
                    center[d],
                    cellMin[d],
                    cellMax[d]);

            const double closestDelta =
                closest - center[d];

            minDistanceSquared +=
                closestDelta *
                closestDelta;

            const double farthestDelta =
                std::max(
                    std::abs(
                        cellMin[d] -
                        center[d]),
                    std::abs(
                        cellMax[d] -
                        center[d]));

            maxDistanceSquared +=
                farthestDelta *
                farthestDelta;
        }

        if(minDistanceSquared <
               radiusSquared &&
           maxDistanceSquared >
               radiusSquared)
        {
            return
                SurfaceRelation::Cross;
        }

        if(minDistanceSquared ==
               radiusSquared ||
           maxDistanceSquared ==
               radiusSquared)
        {
            return
                SurfaceRelation::Touch;
        }

        return
            SurfaceRelation::None;
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
                "Open box must contain the Sphere bounds.");
        }
    }


    [[nodiscard]]
    BoundingBox readOpenBox() const
    {
        std::cout
            << "Sphere bounds:\n"
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
