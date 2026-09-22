#pragma once

#include "bounding_box.hpp"
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
};


//=============================================================================
// STL geometry construction
//=============================================================================

[[nodiscard]]
STLGeometry constructSTLGeometry(
    stl::FacetTopology topology);

} // namespace ntic::lbm::geometry