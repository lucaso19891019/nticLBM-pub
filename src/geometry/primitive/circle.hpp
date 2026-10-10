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

struct Circle
{
    Point2D center;
    double radius;

    BoundingBox2D bounds;
    BoundingBox2D openBox;

    FlowType flowType;


    Circle(
        const Point2D& inputCenter,
        const double inputRadius,
        const FlowType inputFlowType = FlowType::Internal,
        const BoundingBox2D* inputOpenBox = nullptr)
        : center(inputCenter),
          radius(inputRadius),
          bounds{
              {
                  inputCenter[0] - inputRadius,
                  inputCenter[1] - inputRadius
              },
              {
                  inputCenter[0] + inputRadius,
                  inputCenter[1] + inputRadius
              }
          },
          flowType(inputFlowType)
    {
        if(radius <= 0.0)
        {
            throw std::invalid_argument(
                "Circle radius must be positive.");
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
    BoundingBox2D boundingBox() const noexcept
    {
        return bounds;
    }


    void translate(
        const Point2D* targetPoint = nullptr)
    {
        Point2D displacement{};

        if(flowType == FlowType::Internal)
        {
            if(targetPoint == nullptr)
            {
                throw std::invalid_argument(
                    "Internal Circle translation requires a target point.");
            }

            displacement = {
                (*targetPoint)[0] - bounds.min[0],
                (*targetPoint)[1] - bounds.min[1]
            };
        }
        else if(flowType == FlowType::External)
        {
            if(targetPoint != nullptr)
            {
                throw std::invalid_argument(
                    "External Circle translation does not accept a target point.");
            }

            displacement = {
                -openBox.min[0],
                -openBox.min[1]
            };
        }
        else
        {
            throw std::invalid_argument(
                "Unsupported Circle flow type.");
        }

        center[0] += displacement[0];
        center[1] += displacement[1];

        bounds.min[0] += displacement[0];
        bounds.min[1] += displacement[1];

        bounds.max[0] += displacement[0];
        bounds.max[1] += displacement[1];

        if(flowType == FlowType::External)
        {
            openBox.min[0] += displacement[0];
            openBox.min[1] += displacement[1];

            openBox.max[0] += displacement[0];
            openBox.max[1] += displacement[1];
        }
    }



    template <lattice::LatticeType Type>
    void analysis(
        GeometryAnalysis2D<Type>& analysis) const
    {
        interiorAreaAnalysis(analysis);
        boundaryQAnalysis(analysis);
    }



    template <lattice::LatticeType Type>
    void boundaryQAnalysis(
        GeometryAnalysis2D<Type>& analysis) const
    {
        using Lattice =
            typename GeometryAnalysis2D<Type>::Lattice;
    
        constexpr std::size_t nLinks =
            Lattice::nStencils - 1;
    
        const std::size_t boundaryCount =
            analysis.boundaryX.size();
    
        if(analysis.boundaryY.size() != boundaryCount)
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
    
        // AoS: each Boundary Cell stores Q-1 values.
        // -1 means no intersection on the lattice link.
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
    
            const double C =
                px * px + py * py - radius * radius;
    
            for(std::size_t i = 1;
                i < Lattice::nStencils;
                ++i)
            {
                const double dx =
                    h * static_cast<double>(Lattice::ex[i]);
    
                const double dy =
                    h * static_cast<double>(Lattice::ey[i]);
    
                const double A =
                    dx * dx + dy * dy;
    
                const double B =
                    2.0 * (px * dx + py * dy);
    
                const double discriminant =
                    B * B - 4.0 * A * C;
    
                if(discriminant < 0.0)
                {
                    continue;
                }
    
                const double sqrtD =
                    std::sqrt(discriminant);
    
                // Numerically stable quadratic roots.
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
    
                analysis.q[
                    b * nLinks + (i - 1)] = qValue;
            }
        }
    }



    template <lattice::LatticeType Type>
    void interiorAreaAnalysis(
        GeometryAnalysis2D<Type>& analysis) const
    {
        const double gridSpacing = analysis.gridSpacing;

        if(gridSpacing <= 0.0)
        {
            throw std::invalid_argument(
                "Grid spacing must be positive.");
        }

        analysis.domain =
            flowType == FlowType::Internal
                ? bounds
                : openBox;

        const BoundingBox2D& domain = analysis.domain;

        analysis.nx = static_cast<std::size_t>(
            std::ceil(domain.width() / gridSpacing));

        analysis.ny = static_cast<std::size_t>(
            std::ceil(domain.height() / gridSpacing));

        const std::size_t nx = analysis.nx;
        const std::size_t ny = analysis.ny;
        const std::size_t cellCount = nx * ny;

        std::vector<CellType> cellTypes(
            cellCount, CellType::Dry);

        #pragma omp parallel for schedule(static)
        for(std::ptrdiff_t index = 0;
            index < static_cast<std::ptrdiff_t>(cellCount);
            ++index)
        {
            const std::size_t cellID =
                static_cast<std::size_t>(index);

            const std::size_t i = cellID % nx;
            const std::size_t j = cellID / nx;

            const Point2D point
            {
                domain.min[0] +
                    (static_cast<double>(i) + 0.5) *
                    gridSpacing,

                domain.min[1] +
                    (static_cast<double>(j) + 0.5) *
                    gridSpacing
            };

            const bool inside = pointInside(point);

            const bool wet =
                flowType == FlowType::Internal
                    ? inside
                    : !inside;

            if(!wet)
            {
                cellTypes[cellID] = CellType::Dry;
                continue;
            }

            if(distanceToBoundary(point) <= 0.0)
            {
                cellTypes[cellID] = CellType::Dry;
                continue;
            }

            const Point2D stencilMin
            {
                point[0] - gridSpacing,
                point[1] - gridSpacing
            };

            const Point2D stencilMax
            {
                point[0] + gridSpacing,
                point[1] + gridSpacing
            };

            const SurfaceRelation relation =
                surfaceRelation(stencilMin, stencilMax);

            cellTypes[cellID] =
                relation == SurfaceRelation::None
                    ? CellType::Interior
                    : CellType::Boundary;
        }

        auto& boundaryX = analysis.boundaryX;
        auto& boundaryY = analysis.boundaryY;

        boundaryX.resize(cellCount);
        boundaryY.resize(cellCount);

        std::size_t boundaryCount = 0;

        #pragma omp parallel for schedule(static)
        for(std::ptrdiff_t index = 0;
            index < static_cast<std::ptrdiff_t>(cellCount);
            ++index)
        {
            const std::size_t cellID =
                static_cast<std::size_t>(index);

            if(cellTypes[cellID] != CellType::Boundary)
            {
                continue;
            }

            const std::size_t i = cellID % nx;
            const std::size_t j = cellID / nx;

            std::size_t slot = 0;

            #pragma omp atomic capture
            slot = boundaryCount++;

            boundaryX[slot] =
                domain.min[0] +
                (static_cast<double>(i) + 0.5) *
                gridSpacing;

            boundaryY[slot] =
                domain.min[1] +
                (static_cast<double>(j) + 0.5) *
                gridSpacing;
        }

        boundaryX.resize(boundaryCount);
        boundaryY.resize(boundaryCount);

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
        const Point2D& point) const noexcept
    {
        const double dx =
            point[0] - center[0];

        const double dy =
            point[1] - center[1];

        return
            dx * dx +
            dy * dy <
            radius * radius;
    }


    [[nodiscard]]
    double distanceToBoundary(
        const Point2D& point) const noexcept
    {
        const double dx =
            point[0] - center[0];

        const double dy =
            point[1] - center[1];

        const double distance =
            std::sqrt(
                dx * dx +
                dy * dy);

        return
            std::abs(
                distance - radius);
    }


    [[nodiscard]]
    SurfaceRelation surfaceRelation(
        const Point2D& cellMin,
        const Point2D& cellMax) const noexcept
    {
        const double closestX =
            std::clamp(
                center[0],
                cellMin[0],
                cellMax[0]);

        const double closestY =
            std::clamp(
                center[1],
                cellMin[1],
                cellMax[1]);

        const double closestDX =
            closestX - center[0];

        const double closestDY =
            closestY - center[1];

        const double minDistanceSquared =
            closestDX * closestDX +
            closestDY * closestDY;

        const double farthestDX =
            std::max(
                std::abs(
                    cellMin[0] - center[0]),
                std::abs(
                    cellMax[0] - center[0]));

        const double farthestDY =
            std::max(
                std::abs(
                    cellMin[1] - center[1]),
                std::abs(
                    cellMax[1] - center[1]));

        const double maxDistanceSquared =
            farthestDX * farthestDX +
            farthestDY * farthestDY;

        const double radiusSquared =
            radius * radius;

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
        const BoundingBox2D& box) const
    {
        if(box.width() <= 0.0 ||
           box.height() <= 0.0)
        {
            throw std::invalid_argument(
                "Open box dimensions must be positive.");
        }

        if(!box.contains(bounds))
        {
            throw std::invalid_argument(
                "Open box must contain the Circle bounds.");
        }
    }


    [[nodiscard]]
    BoundingBox2D readOpenBox() const
    {
        std::cout
            << "Circle bounds:\n"
            << "  min = ("
            << bounds.min[0] << ", "
            << bounds.min[1] << ")\n"
            << "  max = ("
            << bounds.max[0] << ", "
            << bounds.max[1] << ")\n";

        Point2D origin{};
        Point2D size{};

        std::cout
            << "Enter open box origin x y: ";

        if(!(std::cin >>
             origin[0] >>
             origin[1]))
        {
            throw std::runtime_error(
                "Failed to read open box origin.");
        }

        std::cout
            << "Enter open box size x y: ";

        if(!(std::cin >>
             size[0] >>
             size[1]))
        {
            throw std::runtime_error(
                "Failed to read open box size.");
        }

        const BoundingBox2D box{
            origin,
            {
                origin[0] + size[0],
                origin[1] + size[1]
            }
        };

        validateOpenBox(box);

        return box;
    }
};

} // namespace ntic::lbm::geometry
