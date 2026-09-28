#pragma once

#include "axis.hpp"
#include "bounding_box.hpp"
#include "cell_type.hpp"
#include "flow_type.hpp"
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


    void interiorAreaAnalysis(
        const double gridSpacing,
        std::vector<CellType>& cellTypes,
        std::vector<double>& boundaryX,
        std::vector<double>& boundaryY,
        std::vector<double>& boundaryZ) const
    {
        if(gridSpacing <= 0.0)
        {
            throw std::invalid_argument(
                "Grid spacing must be positive.");
        }

        const BoundingBox& domain =
            flowType == FlowType::Internal
                ? bounds
                : openBox;

        const std::size_t nx =
            static_cast<std::size_t>(
                std::ceil(
                    domain.width() /
                    gridSpacing));

        const std::size_t ny =
            static_cast<std::size_t>(
                std::ceil(
                    domain.height() /
                    gridSpacing));

        const std::size_t nz =
            static_cast<std::size_t>(
                std::ceil(
                    domain.depth() /
                    gridSpacing));

        const std::size_t cellCount =
            nx * ny * nz;

        cellTypes.assign(
            cellCount,
            CellType::Dry);

        const double halfSpacing =
            0.5 * gridSpacing;

        const double halfDiagonal =
            std::sqrt(3.0) *
            halfSpacing;

        #pragma omp parallel for schedule(static)
        for(std::size_t index = 0;
            index < cellCount;
            ++index)
        {
            const std::size_t i =
                index % nx;

            const std::size_t tmp =
                index / nx;

            const std::size_t j =
                tmp % ny;

            const std::size_t k =
                tmp / ny;

            const Point cellCenter{
                domain.min[0] +
                    (static_cast<double>(i) + 0.5) *
                    gridSpacing,

                domain.min[1] +
                    (static_cast<double>(j) + 0.5) *
                    gridSpacing,

                domain.min[2] +
                    (static_cast<double>(k) + 0.5) *
                    gridSpacing
            };

            const bool centerInside =
                pointInside(
                    cellCenter);

            const bool centerFluid =
                flowType == FlowType::Internal
                    ? centerInside
                    : !centerInside;

            if(distanceToBoundary(cellCenter) >
               halfDiagonal)
            {
                cellTypes[index] =
                    centerFluid
                        ? CellType::Interior
                        : CellType::Dry;

                continue;
            }

            const Point cellMin{
                cellCenter[0] - halfSpacing,
                cellCenter[1] - halfSpacing,
                cellCenter[2] - halfSpacing
            };

            const Point cellMax{
                cellCenter[0] + halfSpacing,
                cellCenter[1] + halfSpacing,
                cellCenter[2] + halfSpacing
            };

            const SurfaceRelation relation =
                surfaceRelation(
                    cellMin,
                    cellMax);

            if(relation ==
               SurfaceRelation::Cross)
            {
                cellTypes[index] =
                    CellType::Boundary;

                continue;
            }

            if(relation ==
               SurfaceRelation::Touch)
            {
                cellTypes[index] =
                    cellTouchesFluidSide(
                        cellMin,
                        cellMax)
                        ? CellType::Boundary
                        : CellType::Dry;

                continue;
            }

            cellTypes[index] =
                centerFluid
                    ? CellType::Interior
                    : CellType::Dry;
        }

        boundaryX.clear();
        boundaryY.clear();
        boundaryZ.clear();

        for(std::size_t index = 0;
            index < cellCount;
            ++index)
        {
            if(cellTypes[index] !=
               CellType::Boundary)
            {
                continue;
            }

            const std::size_t i =
                index % nx;

            const std::size_t tmp =
                index / nx;

            const std::size_t j =
                tmp % ny;

            const std::size_t k =
                tmp / ny;

            boundaryX.push_back(
                domain.min[0] +
                (static_cast<double>(i) + 0.5) *
                gridSpacing);

            boundaryY.push_back(
                domain.min[1] +
                (static_cast<double>(j) + 0.5) *
                gridSpacing);

            boundaryZ.push_back(
                domain.min[2] +
                (static_cast<double>(k) + 0.5) *
                gridSpacing);
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


    [[nodiscard]]
    bool cellTouchesFluidSide(
        const Point& cellMin,
        const Point& cellMax) const noexcept
    {
        for(int ix = 0;
            ix < 2;
            ++ix)
        {
            for(int iy = 0;
                iy < 2;
                ++iy)
            {
                for(int iz = 0;
                    iz < 2;
                    ++iz)
                {
                    const Point corner{
                        ix == 0
                            ? cellMin[0]
                            : cellMax[0],

                        iy == 0
                            ? cellMin[1]
                            : cellMax[1],

                        iz == 0
                            ? cellMin[2]
                            : cellMax[2]
                    };

                    const bool inside =
                        pointInside(corner);

                    const bool fluid =
                        flowType ==
                        FlowType::Internal
                            ? inside
                            : !inside;

                    if(fluid)
                    {
                        return true;
                    }
                }
            }
        }

        return false;
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