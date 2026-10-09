#pragma once

#include "axis.hpp"
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

struct Cylinder
{
    Point center;
    double radius;
    double length;
    Axis axis;

    BoundingBox bounds;
    BoundingBox openBox;

    FlowType flowType;


    Cylinder(
        const Point& inputCenter,
        const double inputRadius,
        const double inputLength,
        const Axis inputAxis,
        const FlowType inputFlowType = FlowType::Internal,
        const BoundingBox* inputOpenBox = nullptr)
        : center(inputCenter),
          radius(inputRadius),
          length(inputLength),
          axis(inputAxis),
          flowType(inputFlowType)
    {
        if(radius <= 0.0)
        {
            throw std::invalid_argument(
                "Cylinder radius must be positive.");
        }

        if(length <= 0.0)
        {
            throw std::invalid_argument(
                "Cylinder length must be positive.");
        }

        bounds =
            makeBounds();

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
                    "Internal Cylinder translation requires a target point.");
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
                    "External Cylinder translation does not accept a target point.");
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
                "Unsupported Cylinder flow type.");
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
    void boundaryQAnalysis(GeometryAnalysis3D<Type>&) const
    {
        // Reserved for lattice-link q analysis.
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

    enum class SurfaceRelation
    {
        None,
        Touch,
        Cross
    };


    [[nodiscard]]
    BoundingBox makeBounds() const noexcept
    {
        const double halfLength =
            0.5 * length;

        switch(axis)
        {
            case Axis::X:
                return {
                    {
                        center[0] - halfLength,
                        center[1] - radius,
                        center[2] - radius
                    },
                    {
                        center[0] + halfLength,
                        center[1] + radius,
                        center[2] + radius
                    }
                };

            case Axis::Y:
                return {
                    {
                        center[0] - radius,
                        center[1] - halfLength,
                        center[2] - radius
                    },
                    {
                        center[0] + radius,
                        center[1] + halfLength,
                        center[2] + radius
                    }
                };

            case Axis::Z:
                return {
                    {
                        center[0] - radius,
                        center[1] - radius,
                        center[2] - halfLength
                    },
                    {
                        center[0] + radius,
                        center[1] + radius,
                        center[2] + halfLength
                    }
                };
        }

        return {};
    }


    void coordinates(
        const Point& point,
        double& axial,
        double& radial0,
        double& radial1) const noexcept
    {
        switch(axis)
        {
            case Axis::X:
                axial =
                    point[0] - center[0];

                radial0 =
                    point[1] - center[1];

                radial1 =
                    point[2] - center[2];

                return;

            case Axis::Y:
                axial =
                    point[1] - center[1];

                radial0 =
                    point[0] - center[0];

                radial1 =
                    point[2] - center[2];

                return;

            case Axis::Z:
                axial =
                    point[2] - center[2];

                radial0 =
                    point[0] - center[0];

                radial1 =
                    point[1] - center[1];

                return;
        }
    }


    [[nodiscard]]
    bool pointInside(
        const Point& point) const noexcept
    {
        double axial = 0.0;
        double radial0 = 0.0;
        double radial1 = 0.0;

        coordinates(
            point,
            axial,
            radial0,
            radial1);

        return
            std::abs(axial) <
                0.5 * length &&
            radial0 * radial0 +
            radial1 * radial1 <
                radius * radius;
    }


    [[nodiscard]]
    double distanceToBoundary(
        const Point& point) const noexcept
    {
        double axial = 0.0;
        double radial0 = 0.0;
        double radial1 = 0.0;

        coordinates(
            point,
            axial,
            radial0,
            radial1);

        const double radialDistance =
            std::sqrt(
                radial0 * radial0 +
                radial1 * radial1);

        const double qRadial =
            radialDistance - radius;

        const double qAxial =
            std::abs(axial) -
            0.5 * length;

        const double outsideRadial =
            std::max(
                qRadial,
                0.0);

        const double outsideAxial =
            std::max(
                qAxial,
                0.0);

        const double outsideDistance =
            std::sqrt(
                outsideRadial *
                    outsideRadial +
                outsideAxial *
                    outsideAxial);

        const double insideDistance =
            std::min(
                std::max(
                    qRadial,
                    qAxial),
                0.0);

        return
            std::abs(
                outsideDistance +
                insideDistance);
    }


    [[nodiscard]]
    SurfaceRelation surfaceRelation(
        const Point& cellMin,
        const Point& cellMax) const noexcept
    {
        std::size_t axialDimension = 2;
        std::size_t radialDimension0 = 0;
        std::size_t radialDimension1 = 1;

        if(axis == Axis::X)
        {
            axialDimension = 0;
            radialDimension0 = 1;
            radialDimension1 = 2;
        }
        else if(axis == Axis::Y)
        {
            axialDimension = 1;
            radialDimension0 = 0;
            radialDimension1 = 2;
        }

        const double axialMin =
            center[axialDimension] -
            0.5 * length;

        const double axialMax =
            center[axialDimension] +
            0.5 * length;

        const bool axialClosedOverlap =
            cellMax[axialDimension] >=
                axialMin &&
            cellMin[axialDimension] <=
                axialMax;

        if(!axialClosedOverlap)
        {
            return
                SurfaceRelation::None;
        }

        const double closest0 =
            std::clamp(
                center[radialDimension0],
                cellMin[radialDimension0],
                cellMax[radialDimension0]);

        const double closest1 =
            std::clamp(
                center[radialDimension1],
                cellMin[radialDimension1],
                cellMax[radialDimension1]);

        const double closestDelta0 =
            closest0 -
            center[radialDimension0];

        const double closestDelta1 =
            closest1 -
            center[radialDimension1];

        const double minRadialSquared =
            closestDelta0 *
                closestDelta0 +
            closestDelta1 *
                closestDelta1;

        const double farthestDelta0 =
            std::max(
                std::abs(
                    cellMin[radialDimension0] -
                    center[radialDimension0]),
                std::abs(
                    cellMax[radialDimension0] -
                    center[radialDimension0]));

        const double farthestDelta1 =
            std::max(
                std::abs(
                    cellMin[radialDimension1] -
                    center[radialDimension1]),
                std::abs(
                    cellMax[radialDimension1] -
                    center[radialDimension1]));

        const double maxRadialSquared =
            farthestDelta0 *
                farthestDelta0 +
            farthestDelta1 *
                farthestDelta1;

        const double radiusSquared =
            radius * radius;

        const bool radialClosedOverlap =
            minRadialSquared <=
            radiusSquared;

        if(!radialClosedOverlap)
        {
            return
                SurfaceRelation::None;
        }

        const bool axialOpenOverlap =
            cellMax[axialDimension] >
                axialMin &&
            cellMin[axialDimension] <
                axialMax;

        const bool radialOpenOverlap =
            minRadialSquared <
                radiusSquared;

        const bool cellInsideClosed =
            cellMin[axialDimension] >=
                axialMin &&
            cellMax[axialDimension] <=
                axialMax &&
            maxRadialSquared <=
                radiusSquared;

        const bool cellInsideOpen =
            cellMin[axialDimension] >
                axialMin &&
            cellMax[axialDimension] <
                axialMax &&
            maxRadialSquared <
                radiusSquared;

        if(cellInsideOpen)
        {
            return
                SurfaceRelation::None;
        }

        if(axialOpenOverlap &&
           radialOpenOverlap &&
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
                "Open box must contain the Cylinder bounds.");
        }
    }


    [[nodiscard]]
    BoundingBox readOpenBox() const
    {
        std::cout
            << "Cylinder bounds:\n"
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
