#include "stl_geometry.hpp"

#include <cstddef>
#include <stdexcept>
#include <utility>


namespace ntic::lbm::geometry
{

namespace
{

//=============================================================================
// STL bounding box
//=============================================================================

BoundingBox makeSTLBounds(
    const stl::FacetTopology& topology)
{
    const auto& vertices =
        topology.geometry.vertices;


    const std::size_t vertexCount =
        vertices.size();


    BoundingBox bounds;


#pragma omp parallel
    {
        BoundingBox localBounds;

        bool hasVertex =
            false;


#pragma omp for schedule(static)
        for(std::ptrdiff_t index = 0;
            index <
                static_cast<std::ptrdiff_t>(
                    vertexCount);
            ++index)
        {
            localBounds.expand(
                vertices[
                    static_cast<std::size_t>(
                        index)]);

            hasVertex =
                true;
        }


        if(hasVertex)
        {
#pragma omp critical
            {
                bounds.expand(
                    localBounds.min);

                bounds.expand(
                    localBounds.max);
            }
        }
    }


    return bounds;
}

} // namespace


//=============================================================================
// STL geometry construction
//=============================================================================

STLGeometry::STLGeometry(
    stl::FacetTopology inputTopology,
    const FlowType flowType)
    :
    topology(
        std::move(inputTopology))
{
    if(topology.geometry.vertices.empty())
    {
        throw std::runtime_error(
            "Cannot construct empty STL geometry.");
    }


    bounds =
        makeSTLBounds(
            topology);


    analyzeContainment();


    interpretFlow(
        flowType);
}

} // namespace ntic::lbm::geometry