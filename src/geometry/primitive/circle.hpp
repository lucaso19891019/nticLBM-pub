#pragma once

#include "point.hpp"
#include "bounding_box.hpp"

#include <stdexcept>

namespace ntic::lbm::geometry
{

//=============================================================================
// Circle
//=============================================================================
//
// Two-dimensional circular primitive geometry.
//
// The circle is defined by:
//
//   - center
//   - radius
//
// The circle lies in the z = 0 plane.
//
// Boundary points are considered inside.
//
//=============================================================================

struct Circle
{
    Point center;
    double radius;


    bool contains(const Point& point) const
    {
        if(point[2] != 0.0)
        {
            throw std::invalid_argument(
                "Circle contains() requires point z = 0.");
        }


        const double dx =
            point[0] - center[0];

        const double dy =
            point[1] - center[1];


        return
            dx * dx +
            dy * dy
            <= radius * radius;
    }


    BoundingBox boundingBox() const noexcept
    {
        return {
            {
                center[0] - radius,
                center[1] - radius,
                0.0
            },
            {
                center[0] + radius,
                center[1] + radius,
                0.0
            }
        };
    }
};

} // namespace ntic::lbm::geometry