#pragma once

#include "flow_type.hpp"
#include "stl_geometry.hpp"


namespace ntic::lbm::geometry
{

//=============================================================================
// STL flow interpretation
//=============================================================================

void interpretSTLFlow(
    STLGeometry& geometry,
    FlowType flowType);

} // namespace ntic::lbm::geometry