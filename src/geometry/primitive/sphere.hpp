#pragma once

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

        const double radiusSquared =
            radius * radius;

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
                    cellMax,
                    radiusSquared);

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
        const Point& cellMax,
        const double radiusSquared) const noexcept
    {
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