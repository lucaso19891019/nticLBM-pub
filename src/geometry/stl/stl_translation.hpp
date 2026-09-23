#pragma once

#include "bounding_box.hpp"
#include "point.hpp"
#include "stl_geometry.hpp"


namespace ntic::lbm::geometry
{

//=============================================================================
// STL translation
//=============================================================================

void translate(
    STLGeometry& geometry,
    const Point& displacement);


void translateInternal(
    STLGeometry& geometry,
    const Point& targetPoint);


void translateExternal(
    STLGeometry& geometry,
    BoundingBox& openBox);

} // namespace ntic::lbm::geometry