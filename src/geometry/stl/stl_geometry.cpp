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
    const STLGeometry& geometry)
{
    const std::size_t vertexCount =
        geometry.topology.geometry.vertices.size();


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
                geometry.topology.geometry.vertices[
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

STLGeometry constructSTLGeometry(
    stl::FacetTopology topology)
{
    if(topology.geometry.vertices.empty())
    {
        throw std::runtime_error(
            "Cannot construct empty STL geometry.");
    }


    STLGeometry geometry;

    geometry.topology =
        std::move(topology);


    geometry.bounds =
        makeSTLBounds(
            geometry);


    return geometry;
}

} // namespace ntic::lbm::geometry