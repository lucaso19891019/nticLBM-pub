#include "stl_geometry.hpp"

#include <cstddef>
#include <stdexcept>


namespace ntic::lbm::geometry
{

namespace
{

//=============================================================================
// Translate point
//=============================================================================

void translatePoint(
    Point& point,
    const Point& displacement)
{
    point[0] +=
        displacement[0];

    point[1] +=
        displacement[1];

    point[2] +=
        displacement[2];
}


//=============================================================================
// Translate bounding box
//=============================================================================

void translateBoundingBox(
    BoundingBox& bounds,
    const Point& displacement)
{
    translatePoint(
        bounds.min,
        displacement);

    translatePoint(
        bounds.max,
        displacement);
}


//=============================================================================
// Translate canonical vertices
//=============================================================================

void translateVertices(
    STLGeometry& geometry,
    const Point& displacement)
{
    auto& vertices =
        geometry.topology.geometry.vertices;


    const std::size_t vertexCount =
        vertices.size();


#pragma omp parallel for schedule(static)
    for(std::ptrdiff_t index = 0;
        index <
            static_cast<std::ptrdiff_t>(
                vertexCount);
        ++index)
    {
        translatePoint(
            vertices[
                static_cast<std::size_t>(
                    index)],
            displacement);
    }
}


//=============================================================================
// Translate facet centroids
//=============================================================================

void translateFacetCentroids(
    STLGeometry& geometry,
    const Point& displacement)
{
    auto& facetGeometry =
        geometry.topology.facetGeometry;


    const std::size_t facetCount =
        facetGeometry.size();


#pragma omp parallel for schedule(static)
    for(std::ptrdiff_t index = 0;
        index <
            static_cast<std::ptrdiff_t>(
                facetCount);
        ++index)
    {
        translatePoint(
            facetGeometry[
                static_cast<std::size_t>(
                    index)].centroid,
            displacement);
    }
}


//=============================================================================
// Translate component bounds
//=============================================================================

void translateComponentBounds(
    STLGeometry& geometry,
    const Point& displacement)
{
    auto& components =
        geometry.topology.components;


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


        bounds.min[0] +=
            displacement[0];

        bounds.min[1] +=
            displacement[1];

        bounds.min[2] +=
            displacement[2];


        bounds.max[0] +=
            displacement[0];

        bounds.max[1] +=
            displacement[1];

        bounds.max[2] +=
            displacement[2];


        bounds.center[0] +=
            displacement[0];

        bounds.center[1] +=
            displacement[1];

        bounds.center[2] +=
            displacement[2];
    }
}


//=============================================================================
// Translate STL geometry
//=============================================================================

void translateGeometry(
    STLGeometry& geometry,
    const Point& displacement)
{
    translateVertices(
        geometry,
        displacement);


    translateFacetCentroids(
        geometry,
        displacement);


    translateComponentBounds(
        geometry,
        displacement);


    translateBoundingBox(
        geometry.bounds,
        displacement);
}

} // namespace


//=============================================================================
// STL translation
//=============================================================================

void STLGeometry::translate(
    const Point* targetPoint,
    BoundingBox* openBox)
{
    if(flowType ==
       FlowType::Internal)
    {
        if(targetPoint ==
           nullptr)
        {
            throw std::runtime_error(
                "Internal STL translation requires a target point.");
        }


        const Point displacement =
        {{
            (*targetPoint)[0] -
                bounds.min[0],

            (*targetPoint)[1] -
                bounds.min[1],

            (*targetPoint)[2] -
                bounds.min[2]
        }};


        translateGeometry(
            *this,
            displacement);


        return;
    }


    if(flowType ==
       FlowType::External)
    {
        if(openBox ==
           nullptr)
        {
            throw std::runtime_error(
                "External STL translation requires an open box.");
        }


        const Point displacement =
        {{
            -openBox->min[0],
            -openBox->min[1],
            -openBox->min[2]
        }};


        translateGeometry(
            *this,
            displacement);


        translateBoundingBox(
            *openBox,
            displacement);


        return;
    }


    throw std::runtime_error(
        "Unsupported STL flow type.");
}

} // namespace ntic::lbm::geometry