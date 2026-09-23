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

} // namespace


//=============================================================================
// STL translation
//=============================================================================

void STLGeometry::translate(
    const Point* targetPoint,
    BoundingBox* openBox)
{
    Point displacement;


    if(flowType ==
       FlowType::Internal)
    {
        if(targetPoint ==
           nullptr)
        {
            throw std::runtime_error(
                "Internal STL translation requires a target point.");
        }


        displacement =
        {{
            (*targetPoint)[0] -
                bounds.min[0],

            (*targetPoint)[1] -
                bounds.min[1],

            (*targetPoint)[2] -
                bounds.min[2]
        }};
    }
    else if(flowType ==
            FlowType::External)
    {
        if(openBox ==
           nullptr)
        {
            throw std::runtime_error(
                "External STL translation requires an open box.");
        }


        displacement =
        {{
            -openBox->min[0],
            -openBox->min[1],
            -openBox->min[2]
        }};
    }
    else
    {
        throw std::runtime_error(
            "Unsupported STL flow type.");
    }


    //-------------------------------------------------------------------------
    // Canonical vertices
    //-------------------------------------------------------------------------

    auto& vertices =
        topology.geometry.vertices;


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


    //-------------------------------------------------------------------------
    // Facet centroids
    //-------------------------------------------------------------------------

    auto& facetGeometry =
        topology.facetGeometry;


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


    //-------------------------------------------------------------------------
    // Component bounds
    //-------------------------------------------------------------------------

    auto& components =
        topology.components;


    const std::size_t componentCount =
        components.size();


#pragma omp parallel for schedule(static)
    for(std::ptrdiff_t index = 0;
        index <
            static_cast<std::ptrdiff_t>(
                componentCount);
        ++index)
    {
        auto& componentBounds =
            components[
                static_cast<std::size_t>(
                    index)].bounds;


        componentBounds.min[0] +=
            displacement[0];

        componentBounds.min[1] +=
            displacement[1];

        componentBounds.min[2] +=
            displacement[2];


        componentBounds.max[0] +=
            displacement[0];

        componentBounds.max[1] +=
            displacement[1];

        componentBounds.max[2] +=
            displacement[2];


        componentBounds.center[0] +=
            displacement[0];

        componentBounds.center[1] +=
            displacement[1];

        componentBounds.center[2] +=
            displacement[2];
    }


    //-------------------------------------------------------------------------
    // Global bounds
    //-------------------------------------------------------------------------

    translateBoundingBox(
        bounds,
        displacement);


    //-------------------------------------------------------------------------
    // External open box
    //-------------------------------------------------------------------------

    if(flowType ==
       FlowType::External)
    {
        translateBoundingBox(
            *openBox,
            displacement);
    }
}

} // namespace ntic::lbm::geometry