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

            const Point center{
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
                pointInside(center);

            const bool centerFluid =
                flowType == FlowType::Internal
                    ? centerInside
                    : !centerInside;

            if(distanceToBoundary(center) >
               halfDiagonal)
            {
                cellTypes[index] =
                    centerFluid
                        ? CellType::Interior
                        : CellType::Dry;

                continue;
            }

            const Point cellMin{
                domain.min[0] +
                    static_cast<double>(i) *
                    gridSpacing,

                domain.min[1] +
                    static_cast<double>(j) *
                    gridSpacing,

                domain.min[2] +
                    static_cast<double>(k) *
                    gridSpacing
            };

            const Point cellMax{
                domain.min[0] +
                    (static_cast<double>(i) + 1.0) *
                    gridSpacing,

                domain.min[1] +
                    (static_cast<double>(j) + 1.0) *
                    gridSpacing,

                domain.min[2] +
                    (static_cast<double>(k) + 1.0) *
                    gridSpacing
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