#pragma once

#include "bounding_box.hpp"
#include "flow_type.hpp"
#include "point.hpp"
#include "stl_topology.hpp"

#include <cstddef>
#include <limits>
#include <vector>


namespace ntic::lbm::geometry
{

//=============================================================================
// STL component containment
//=============================================================================

struct STLComponentContainment
{
    std::size_t parent =
        std::numeric_limits<std::size_t>::max();

    std::vector<std::size_t> children;

    std::size_t root =
        std::numeric_limits<std::size_t>::max();

    std::size_t level =
        0;
};


//=============================================================================
// STL component flow
//=============================================================================

enum class FluidSide
{
    Inside,
    Outside
};


struct STLComponentFlow
{
    bool active =
        false;

    FluidSide fluidSide =
        FluidSide::Outside;
};


//=============================================================================
// STLGeometry
//=============================================================================

struct STLGeometry
{
    stl::FacetTopology topology;

    BoundingBox bounds;

    std::vector<STLComponentContainment> containment;

    std::vector<STLComponentFlow> flow;

    FlowType flowType;


    STLGeometry(
        stl::FacetTopology topology,
        FlowType flowType);


    void analyzeContainment();


    void interpretFlow(
        FlowType flowType);


    void translate(
        const Point* targetPoint = nullptr,
        BoundingBox* openBox = nullptr);


    [[nodiscard]]
    bool contains(
        const Point& point) const;
};

} // namespace ntic::lbm::geometry