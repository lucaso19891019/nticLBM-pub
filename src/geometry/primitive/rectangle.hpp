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

struct Rectangle
{
    Point2D min;
    Point2D max;

    BoundingBox2D bounds;
    BoundingBox2D openBox;

    FlowType flowType;


    Rectangle(
        const Point2D& inputMin,
        const Point2D& inputMax,
        const FlowType inputFlowType = FlowType::Internal,
        const BoundingBox2D* inputOpenBox = nullptr)
        : min(inputMin),
          max(inputMax),
          bounds{inputMin, inputMax},
          flowType(inputFlowType)
    {
        if(max[0] <= min[0] ||
           max[1] <= min[1])
        {
            throw std::invalid_argument(
                "Rectangle dimensions must be positive.");
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
                    "Internal Rectangle translation requires a target point.");
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
                    "External Rectangle translation does not accept a target point.");
            }

            displacement = {
                -openBox.min[0],
                -openBox.min[1]
            };
        }
        else
        {
            throw std::invalid_argument(
                "Unsupported Rectangle flow type.");
        }

        min[0] += displacement[0];
        min[1] += displacement[1];

        max[0] += displacement[0];
        max[1] += displacement[1];

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


    void interiorAreaAnalysis(
        const double gridSpacing,
        std::vector<CellType>& cellTypes,
        std::vector<double>& boundaryX,
        std::vector<double>& boundaryY) const
    {
        if(gridSpacing <= 0.0)
        {
            throw std::invalid_argument(
                "Grid spacing must be positive.");
        }

        const BoundingBox2D& domain =
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

        const std::size_t cellCount =
            nx * ny;

        cellTypes.assign(
            cellCount,
            CellType::Dry);

        const double halfSpacing =
            0.5 * gridSpacing;

        const double halfDiagonal =
            std::sqrt(2.0) *
            halfSpacing;

        #pragma omp parallel for schedule(static)
        for(std::size_t index = 0;
            index < cellCount;
            ++index)
        {
            const std::size_t i =
                index % nx;

            const std::size_t j =
                index / nx;

            const Point2D center{
                domain.min[0] +
                    (static_cast<double>(i) + 0.5) *
                    gridSpacing,

                domain.min[1] +
                    (static_cast<double>(j) + 0.5) *
                    gridSpacing
            };

            const bool centerInside =
                pointInside(center);

            const bool centerFluid =
                flowType == FlowType::Internal
                    ? centerInside
                    : !centerInside;

            const double distance =
                distanceToBoundary(center);

            if(distance > halfDiagonal)
            {
                cellTypes[index] =
                    centerFluid
                        ? CellType::Interior
                        : CellType::Dry;

                continue;
            }

            const Point2D cellMin{
                center[0] - halfSpacing,
                center[1] - halfSpacing
            };

            const Point2D cellMax{
                center[0] + halfSpacing,
                center[1] + halfSpacing
            };

            const SurfaceRelation relation =
                surfaceRelation(
                    cellMin,
                    cellMax);

            if(relation == SurfaceRelation::Cross)
            {
                cellTypes[index] =
                    CellType::Boundary;

                continue;
            }

            if(relation == SurfaceRelation::Touch)
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

            const std::size_t j =
                index / nx;

            boundaryX.push_back(
                domain.min[0] +
                (static_cast<double>(i) + 0.5) *
                gridSpacing);

            boundaryY.push_back(
                domain.min[1] +
                (static_cast<double>(j) + 0.5) *
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
        const Point2D& point) const noexcept
    {
        return
            point[0] > min[0] &&
            point[0] < max[0] &&
            point[1] > min[1] &&
            point[1] < max[1];
    }


    [[nodiscard]]
    double distanceToBoundary(
        const Point2D& point) const noexcept
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

        if(dx > 0.0 ||
           dy > 0.0)
        {
            return
                std::sqrt(
                    dx * dx +
                    dy * dy);
        }

        return
            std::min(
                std::min(
                    point[0] - min[0],
                    max[0] - point[0]),
                std::min(
                    point[1] - min[1],
                    max[1] - point[1]));
    }


    [[nodiscard]]
    SurfaceRelation surfaceRelation(
        const Point2D& cellMin,
        const Point2D& cellMax) const noexcept
    {
        const bool closedOverlap =
            cellMax[0] >= min[0] &&
            cellMin[0] <= max[0] &&
            cellMax[1] >= min[1] &&
            cellMin[1] <= max[1];

        if(!closedOverlap)
        {
            return
                SurfaceRelation::None;
        }

        const bool cellInsideOpen =
            cellMin[0] > min[0] &&
            cellMax[0] < max[0] &&
            cellMin[1] > min[1] &&
            cellMax[1] < max[1];

        if(cellInsideOpen)
        {
            return
                SurfaceRelation::None;
        }

        const bool cellInsideClosed =
            cellMin[0] >= min[0] &&
            cellMax[0] <= max[0] &&
            cellMin[1] >= min[1] &&
            cellMax[1] <= max[1];

        const bool openOverlap =
            cellMax[0] > min[0] &&
            cellMin[0] < max[0] &&
            cellMax[1] > min[1] &&
            cellMin[1] < max[1];

        if(openOverlap &&
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
        const Point2D& cellMin,
        const Point2D& cellMax) const noexcept
    {
        const Point2D corners[4]{
            {cellMin[0], cellMin[1]},
            {cellMax[0], cellMin[1]},
            {cellMin[0], cellMax[1]},
            {cellMax[0], cellMax[1]}
        };

        for(const Point2D& corner : corners)
        {
            const bool inside =
                pointInside(corner);

            const bool fluid =
                flowType == FlowType::Internal
                    ? inside
                    : !inside;

            if(fluid)
            {
                return true;
            }
        }

        return false;
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
                "Open box must contain the Rectangle bounds.");
        }
    }


    [[nodiscard]]
    BoundingBox2D readOpenBox() const
    {
        std::cout
            << "Rectangle bounds:\n"
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