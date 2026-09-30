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


    template <typename LatticeModel>
    void analysis(
        GeometryAnalysis2D<LatticeModel>& analysis) const
    {
        interiorAreaAnalysis(
            analysis.gridSpacing,
            analysis.domain,
            analysis.nx,
            analysis.ny,
            analysis.scalar,
            analysis.boundaryX,
            analysis.boundaryY);

        boundaryQAnalysis(
            analysis);
    }


    template <typename LatticeModel>
    void boundaryQAnalysis(
        GeometryAnalysis2D<LatticeModel>&) const
    {
        // Reserved for lattice-link q analysis.
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


    #pragma omp parallel for schedule(static)
        for(std::ptrdiff_t index = 0;
            index <
                static_cast<std::ptrdiff_t>(
                    cellCount);
            ++index)
        {
            const std::size_t cellID =
                static_cast<std::size_t>(
                    index);

            const std::size_t i =
                cellID %
                nx;

            const std::size_t j =
                cellID /
                nx;


            const Point2D point
            {
                domain.min[0] +
                    (static_cast<double>(i) + 0.5) *
                    gridSpacing,

                domain.min[1] +
                    (static_cast<double>(j) + 0.5) *
                    gridSpacing
            };


            const bool inside =
                pointInside(
                    point);

            const bool wet =
                flowType == FlowType::Internal
                    ? inside
                    : !inside;


            if(!wet)
            {
                cellTypes[cellID] =
                    CellType::Dry;

                continue;
            }


            if(distanceToBoundary(point) <=
            std::numeric_limits<double>::epsilon())
            {
                cellTypes[cellID] =
                    CellType::Dry;

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
                surfaceRelation(
                    stencilMin,
                    stencilMax);


            if(relation !=
            SurfaceRelation::None)
            {
                cellTypes[cellID] =
                    CellType::Boundary;
            }
            else
            {
                cellTypes[cellID] =
                    CellType::Interior;
            }
        }


        boundaryX.resize(
            cellCount);

        boundaryY.resize(
            cellCount);


        std::size_t boundaryCount =
            0;


    #pragma omp parallel for schedule(static)
        for(std::ptrdiff_t index = 0;
            index <
                static_cast<std::ptrdiff_t>(
                    cellCount);
            ++index)
        {
            const std::size_t cellID =
                static_cast<std::size_t>(
                    index);


            if(cellTypes[cellID] !=
            CellType::Boundary)
            {
                continue;
            }


            const std::size_t i =
                cellID %
                nx;

            const std::size_t j =
                cellID /
                nx;


            std::size_t slot =
                0;

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


        boundaryX.resize(
            boundaryCount);

        boundaryY.resize(
            boundaryCount);
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