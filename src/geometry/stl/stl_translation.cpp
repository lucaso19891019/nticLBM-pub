#include "stl_geometry.hpp"

#include <cstddef>
#include <stdexcept>


namespace ntic::lbm::geometry
{

namespace
{

void translatePoint(
    Point& point,
    const Point& displacement)
{
    point[0] += displacement[0];
    point[1] += displacement[1];
    point[2] += displacement[2];
}


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


void translateVertices(
    STLGeometry& geometry,
    const Point& displacement)
{
    auto& vertices =
        geometry.topology.geometry.vertices;


#pragma omp parallel for schedule(static)
    for(std::ptrdiff_t index = 0;
        index <
            static_cast<std::ptrdiff_t>(
                vertices.size());
        ++index)
    {
        translatePoint(
            vertices[
                static_cast<std::size_t>(
                    index)],
            displacement);
    }
}


void translateFacetCentroids(
    STLGeometry& geometry,
    const Point& displacement)
{
    auto& facetGeometry =
        geometry.topology.facetGeometry;


#pragma omp parallel for schedule(static)
    for(std::ptrdiff_t index = 0;
        index <
            static_cast<std::ptrdiff_t>(
                facetGeometry.size());
        ++index)
    {
        translatePoint(
            facetGeometry[
                static_cast<std::size_t>(
                    index)].centroid,
            displacement);
    }
}


void translateComponentBounds(
    STLGeometry& geometry,
    const Point& displacement)
{
    auto& components =
        geometry.topology.components;


#pragma omp parallel for schedule(static)
    for(std::ptrdiff_t index = 0;
        index <
            static_cast<std::ptrdiff_t>(
                components.size());
        ++index)
    {
        auto& bounds =
            components[
                static_cast<std::size_t>(
                    index)].bounds;


        bounds.min[0] += displacement[0];
        bounds.min[1] += displacement[1];
        bounds.min[2] += displacement[2];

        bounds.max[0] += displacement[0];
        bounds.max[1] += displacement[1];
        bounds.max[2] += displacement[2];

        bounds.center[0] += displacement[0];
        bounds.center[1] += displacement[1];
        bounds.center[2] += displacement[2];
    }
}


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
    const Point* targetPoint)
{
    if(flowType ==
       FlowType::Internal)
    {
        if(targetPoint == nullptr)
        {
            throw std::runtime_error(
                "Internal STL translation requires a target point.");
        }


        const Point displacement =
        {
            (*targetPoint)[0] - bounds.min[0],
            (*targetPoint)[1] - bounds.min[1],
            (*targetPoint)[2] - bounds.min[2]
        };


        translateGeometry(
            *this,
            displacement);

        return;
    }


    if(flowType ==
       FlowType::External)
    {
        if(targetPoint != nullptr)
        {
            throw std::runtime_error(
                "External STL translation does not accept a target point.");
        }


        const Point displacement =
        {
            -openBox.min[0],
            -openBox.min[1],
            -openBox.min[2]
        };


        translateGeometry(
            *this,
            displacement);

        translateBoundingBox(
            openBox,
            displacement);

        return;
    }


    throw std::runtime_error(
        "Unsupported STL flow type for translation.");
}

} // namespace ntic::lbm::geometry
