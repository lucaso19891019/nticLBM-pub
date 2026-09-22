#pragma once

#include "point.hpp"
#include "bounding_box.hpp"

#include "flow_type.hpp"

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


    bool contains(
        const Point& point,
        const FlowType flowType,
        const BoundingBox* openBox = nullptr) const
    {
        const double dx =
            point[0] - center[0];

        const double dy =
            point[1] - center[1];


        const double distanceSquared =
            dx * dx +
            dy * dy;

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


        const bool insideOpenBox =
            point[0] > openBox->min[0] &&
            point[0] < openBox->max[0] &&

            point[1] > openBox->min[1] &&
            point[1] < openBox->max[1];


        return
            insideOpenBox &&
            distanceSquared >
                radiusSquared;
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