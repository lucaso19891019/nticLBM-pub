#pragma once

#include <array>
#include <algorithm>
#include <cstddef>
#include <limits>

#include "point.hpp"

namespace ntic::lbm::geometry
{
//=============================================================================
// BoundingBox
//=============================================================================
//
// Axis-aligned bounding box used by geometry objects.
//
// The bounding box is represented by its minimum and maximum coordinates:
//
//     min = (xmin, ymin, zmin)
//     max = (xmax, ymax, zmax)
//
// A three-dimensional representation is used consistently throughout the
// geometry module.
//
//=============================================================================

struct BoundingBox
{
    Point min{
        std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity()
    };

    Point max{
        -std::numeric_limits<double>::infinity(),
        -std::numeric_limits<double>::infinity(),
        -std::numeric_limits<double>::infinity()
    };


    double width() const noexcept
    {
        return max[0] - min[0];
    }


    double height() const noexcept
    {
        return max[1] - min[1];
    }


    double depth() const noexcept
    {
        return max[2] - min[2];
    }


    bool contains(const Point& point) const noexcept
    {
        return
            point[0] >= min[0] &&
            point[0] <= max[0] &&

            point[1] >= min[1] &&
            point[1] <= max[1] &&

            point[2] >= min[2] &&
            point[2] <= max[2];
    }

    bool contains(
        const BoundingBox& box) const noexcept
    {
        return
            box.min[0] >= min[0] &&
            box.max[0] <= max[0] &&

            box.min[1] >= min[1] &&
            box.max[1] <= max[1] &&

            box.min[2] >= min[2] &&
            box.max[2] <= max[2];
    }

    void expand(const Point& point) noexcept
    {
        for(std::size_t d = 0; d < 3; ++d)
        {
            min[d] =
                std::min(min[d], point[d]);

            max[d] =
                std::max(max[d], point[d]);
        }
    }
};

} // namespace ntic::lbm::geometry