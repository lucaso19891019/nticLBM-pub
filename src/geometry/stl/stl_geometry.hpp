#pragma once

#include "bounding_box.hpp"
#include "flow_type.hpp"

#include "stl_topology.hpp"

namespace ntic::lbm::geometry
{

struct STLGeometry
{
    stl::FacetTopology topology;

    BoundingBox bounds;

    BoundingBox domainBounds;

    FlowType flowType =
        FlowType::Internal;
};


STLGeometry prepareSTLGeometry(
    stl::FacetTopology topology,
    FlowType flowType,
    const BoundingBox& openBox = {},
    const Point& targetSTLMin = {});

} // namespace ntic::lbm::geometry