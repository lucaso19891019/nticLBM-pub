#include "stl_translation.hpp"

#include <cstddef>


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

} // namespace


//=============================================================================
// STL translation
//=============================================================================

void translate(
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


//=============================================================================
// Internal STL translation
//=============================================================================

void translateInternal(
    STLGeometry& geometry,
    const Point& targetPoint)
{
    const Point displacement =
    {
        targetPoint[0] -
            geometry.bounds.min[0],

        targetPoint[1] -
            geometry.bounds.min[1],

        targetPoint[2] -
            geometry.bounds.min[2]
    };


    translate(
        geometry,
        displacement);
}


//=============================================================================
// External STL translation
//=============================================================================

void translateExternal(
    STLGeometry& geometry,
    BoundingBox& openBox)
{
    const Point displacement =
    {
        -openBox.min[0],
        -openBox.min[1],
        -openBox.min[2]
    };


    translate(
        geometry,
        displacement);


    translateBoundingBox(
        openBox,
        displacement);
}

} // namespace ntic::lbm::geometry