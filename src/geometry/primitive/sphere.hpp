#pragma once

#include "point.hpp"
#include "bounding_box.hpp"

namespace ntic::lbm::geometry
{

//=============================================================================
// Sphere
//=============================================================================
//
// Three-dimensional spherical primitive geometry.
//
// The sphere is defined by:
//
//   - center
//   - radius
//
// Boundary points are considered inside.
//
//=============================================================================

struct Sphere
{
    Point center;
    double radius;


    bool contains(const Point& point) const noexcept
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
            dz * dz
            <= radius * radius;
    }


    BoundingBox boundingBox() const noexcept
    {
        return {
            {
                center[0] - radius,
                center[1] - radius,
                center[2] - radius
            },
            {
                center[0] + radius,
                center[1] + radius,
                center[2] + radius
            }
        };
    }
};

} // namespace ntic::lbm::geometry