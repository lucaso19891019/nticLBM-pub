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
    BoundingBox bounds;


    for(const auto& vertex :
        topology.geometry.vertices)
    {
        bounds.expand(vertex);
    }


    return bounds;
}


//=============================================================================
// Translation
//=============================================================================

Point computeTranslation(
    const BoundingBox& bounds,
    const FlowType flowType,
    const Point& targetSTLMin)
{
    Point translation{};


    if(flowType ==
       FlowType::Internal)
    {
        for(std::size_t d = 0;
            d < 3;
            ++d)
        {
            translation[d] =
                -bounds.min[d];
        }
    }
    else
    {
        for(std::size_t d = 0;
            d < 3;
            ++d)
        {
            translation[d] =
                targetSTLMin[d] -
                bounds.min[d];
        }
    }


    return translation;
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
    const std::size_t componentCount =
        components.size();


#pragma omp parallel for schedule(static)
    for(std::ptrdiff_t index = 0;
        index <
            static_cast<std::ptrdiff_t>(
                componentCount);
        ++index)
    {
        auto& bounds =
            components[
                static_cast<std::size_t>(
                    index)].bounds;


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
// Prepare STL geometry
//=============================================================================

STLGeometry prepareSTLGeometry(
    stl::FacetTopology topology,
    const FlowType flowType,
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
        flowType;

    geometry.bounds =
        makeSTLBounds(
            topology);


    const Point translation =
        computeTranslation(
            geometry.bounds,
            flowType,
            targetSTLMin);


    translateTopology(
        topology,
        translation);

    translateBounds(
        geometry.bounds,
        translation);


    if(flowType ==
       FlowType::Internal)
    {
        geometry.domainBounds =
            geometry.bounds;
    }
    else
    {
        if(!openBox.contains(
               geometry.bounds))
        {
            throw std::runtime_error(
                "Translated STL bounding box is not "
                "fully contained in the open box.");
        }


        geometry.domainBounds =
            openBox;
    }


    geometry.topology =
        std::move(topology);


    return geometry;
}

} // namespace ntic::lbm::geometry