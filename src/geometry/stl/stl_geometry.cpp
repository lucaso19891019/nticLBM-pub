#include "stl_geometry.hpp"

#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>


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
    const std::size_t vertexCount =
        topology.geometry.vertices.size();


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
                topology.geometry.vertices[
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


//=============================================================================
// Translate geometric vertices
//=============================================================================

void translateVertices(
    stl::GeometricVertices& geometry,
    const Point& translation)
{
    const std::size_t vertexCount =
        geometry.vertices.size();


#pragma omp parallel for schedule(static)
    for(std::ptrdiff_t index = 0;
        index <
            static_cast<std::ptrdiff_t>(
                vertexCount);
        ++index)
    {
        auto& vertex =
            geometry.vertices[
                static_cast<std::size_t>(
                    index)];


        for(std::size_t d = 0;
            d < 3;
            ++d)
        {
            vertex[d] +=
                translation[d];
        }
    }
}


//=============================================================================
// Translate facet geometry
//=============================================================================

void translateFacetGeometry(
    std::vector<stl::FacetGeometry>& facetGeometry,
    const Point& translation)
{
    const std::size_t facetCount =
        facetGeometry.size();


#pragma omp parallel for schedule(static)
    for(std::ptrdiff_t index = 0;
        index <
            static_cast<std::ptrdiff_t>(
                facetCount);
        ++index)
    {
        auto& centroid =
            facetGeometry[
                static_cast<std::size_t>(
                    index)].centroid;


        for(std::size_t d = 0;
            d < 3;
            ++d)
        {
            centroid[d] +=
                translation[d];
        }
    }
}


//=============================================================================
// Translate component bounds
//=============================================================================

void translateComponents(
    stl::STLComponents& components,
    const Point& translation)
{
    for(auto& component :
        components)
    {
        auto& bounds =
            component.bounds;


        for(std::size_t d = 0;
            d < 3;
            ++d)
        {
            bounds.min[d] +=
                translation[d];

            bounds.max[d] +=
                translation[d];

            bounds.center[d] +=
                translation[d];
        }
    }
}


//=============================================================================
// Translate STL topology
//=============================================================================

void translateTopology(
    stl::FacetTopology& topology,
    const Point& translation)
{
    translateVertices(
        topology.geometry,
        translation);

    translateFacetGeometry(
        topology.facetGeometry,
        translation);

    translateComponents(
        topology.components,
        translation);
}


//=============================================================================
// Translate bounding box
//=============================================================================

void translateBounds(
    BoundingBox& bounds,
    const Point& translation)
{
    for(std::size_t d = 0;
        d < 3;
        ++d)
    {
        bounds.min[d] +=
            translation[d];

        bounds.max[d] +=
            translation[d];
    }
}

} // namespace


//=============================================================================
// Internal STL geometry
//=============================================================================

STLGeometry prepareInternalSTLGeometry(
    stl::FacetTopology topology)
{
    if(topology.geometry.vertices.empty())
    {
        throw std::runtime_error(
            "Cannot prepare empty STL geometry.");
    }


    STLGeometry geometry;

    geometry.flowType =
        FlowType::Internal;

    geometry.bounds =
        makeSTLBounds(
            topology);


    Point translation{};

    for(std::size_t d = 0;
        d < 3;
        ++d)
    {
        translation[d] =
            -geometry.bounds.min[d];
    }


    translateTopology(
        topology,
        translation);

    translateBounds(
        geometry.bounds,
        translation);


    geometry.domainBounds =
        geometry.bounds;

    geometry.topology =
        std::move(topology);


    return geometry;
}


//=============================================================================
// External STL geometry
//=============================================================================

STLGeometry prepareExternalSTLGeometry(
    stl::FacetTopology topology,
    const BoundingBox& openBox,
    const Point& targetSTLMin)
{
    if(topology.geometry.vertices.empty())
    {
        throw std::runtime_error(
            "Cannot prepare empty STL geometry.");
    }


    STLGeometry geometry;

    geometry.flowType =
        FlowType::External;

    geometry.bounds =
        makeSTLBounds(
            topology);


    Point translation{};

    for(std::size_t d = 0;
        d < 3;
        ++d)
    {
        translation[d] =
            targetSTLMin[d] -
            geometry.bounds.min[d];
    }


    translateTopology(
        topology,
        translation);

    translateBounds(
        geometry.bounds,
        translation);


    if(!openBox.contains(
           geometry.bounds))
    {
        throw std::runtime_error(
            "Translated STL bounding box is not "
            "fully contained in the open box.");
    }


    geometry.domainBounds =
        openBox;

    geometry.topology =
        std::move(topology);


    return geometry;
}

} // namespace ntic::lbm::geometry