#pragma once

#include "point.hpp"
#include "axis.hpp"
#include "bounding_box.hpp"

namespace ntic::lbm::geometry
{
//=============================================================================
// Cylinder
//=============================================================================
//
// Axis-aligned three-dimensional cylindrical primitive geometry.
//
// The cylinder is defined by:
//
//   - center
//   - radius
//   - length
//   - axis (X, Y, or Z)
//
// The center is located at the geometric center of the cylinder.
//
// Boundary points are considered inside.
//
//=============================================================================

struct Cylinder
{
    Point center;
    double radius;
    double length;
    Axis axis;


    bool contains(const Point& point) const noexcept
    {
        const double dx =
            point[0] - center[0];

        const double dy =
            point[1] - center[1];

        const double dz =
            point[2] - center[2];

        const double halfLength =
            0.5 * length;


        switch(axis)
        {
            case Axis::X:
                return
                    dx >= -halfLength &&
                    dx <=  halfLength &&
                    dy * dy + dz * dz <= radius * radius;

            case Axis::Y:
                return
                    dy >= -halfLength &&
                    dy <=  halfLength &&
                    dx * dx + dz * dz <= radius * radius;

            case Axis::Z:
                return
                    dz >= -halfLength &&
                    dz <=  halfLength &&
                    dx * dx + dy * dy <= radius * radius;
        }

        return false;
    }


    BoundingBox boundingBox() const noexcept
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
};

} // namespace ntic::lbm::geometry