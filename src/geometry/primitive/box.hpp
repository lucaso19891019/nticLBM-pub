#pragma once

#include "point.hpp"
#include "bounding_box.hpp"

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


    BoundingBox boundingBox() const noexcept
    {
        return {min, max};
    }
};

} // namespace ntic::lbm::geometry