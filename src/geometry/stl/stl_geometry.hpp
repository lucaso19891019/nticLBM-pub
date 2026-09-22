#pragma once

#include "bounding_box.hpp"
#include "flow_type.hpp"
#include "point.hpp"
#include "stl_topology.hpp"


namespace ntic::lbm::geometry
{

//=============================================================================
// STLGeometry
//=============================================================================

struct STLGeometry
{
    stl::FacetTopology topology;

    BoundingBox bounds;

    BoundingBox domainBounds;

    FlowType flowType =
        FlowType::Internal;
};


//=============================================================================
// STL geometry preparation
//=============================================================================

[[nodiscard]]
STLGeometry prepareInternalSTLGeometry(
    stl::FacetTopology topology);


[[nodiscard]]
STLGeometry prepareExternalSTLGeometry(
    stl::FacetTopology topology,
    const BoundingBox& openBox,
    const Point& targetSTLMin);

} // namespace ntic::lbm::geometry