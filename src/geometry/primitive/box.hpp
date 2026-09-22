#pragma once

#include "point.hpp"
#include "bounding_box.hpp"

#include "flow_type.hpp"

#include <stdexcept>

namespace ntic::lbm::geometry
{

//=============================================================================
// Box
//=============================================================================
//
// Axis-aligned three-dimensional box.
//
// The box is defined by its minimum and maximum coordinates.
//
// It provides:
//
//   - point inside/outside classification
//   - axis-aligned bounding box
//
// Boundary points are considered inside.
//
//=============================================================================

struct Box
{
    Point min;
    Point max;


    bool contains(
        const Point& point,
        const FlowType flowType,
        const BoundingBox* openBox = nullptr) const
    {
        const bool strictlyInside =
            point[0] > min[0] &&
            point[0] < max[0] &&

            point[1] > min[1] &&
            point[1] < max[1] &&

            point[2] > min[2] &&
            point[2] < max[2];


        if(flowType ==
        FlowType::Internal)
        {
            return
                strictlyInside;
        }


        if(openBox == nullptr)
        {
            throw std::invalid_argument(
                "External flow requires an open box.");
        }


        const bool onSurface =
            (
                point[0] >= min[0] &&
                point[0] <= max[0] &&

                point[1] >= min[1] &&
                point[1] <= max[1] &&

                point[2] >= min[2] &&
                point[2] <= max[2]
            ) &&
            (
                point[0] == min[0] ||
                point[0] == max[0] ||

                point[1] == min[1] ||
                point[1] == max[1] ||

                point[2] == min[2] ||
                point[2] == max[2]
            );


        const bool strictlyOutside =
            !strictlyInside &&
            !onSurface;


        return
            openBox->strictlyContains(point) &&
            strictlyOutside;
    }


    BoundingBox boundingBox() const noexcept
    {
        return {min, max};
    }
};

} // namespace ntic::lbm::geometry