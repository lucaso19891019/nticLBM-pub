#pragma once

#include "point.hpp"
#include "bounding_box.hpp"

#include "flow_type.hpp"

#include <stdexcept>

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


    bool contains(
        const Point& point,
        const FlowType flowType,
        const BoundingBox* openBox = nullptr) const
    {
        const double dx =
            point[0] - center[0];

        const double dy =
            point[1] - center[1];

        const double dz =
            point[2] - center[2];

        const double distanceSquared =
            dx * dx +
            dy * dy +
            dz * dz;

        const double radiusSquared =
            radius * radius;


        if(flowType ==
        FlowType::Internal)
        {
            return
                distanceSquared <
                radiusSquared;
        }


        if(openBox == nullptr)
        {
            throw std::invalid_argument(
                "External flow requires an open box.");
        }


        return
            openBox->strictlyContains(point) &&
            distanceSquared >
                radiusSquared;
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