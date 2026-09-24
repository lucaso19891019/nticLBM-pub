#include "stl_geometry.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>


namespace ntic::lbm::geometry
{

namespace
{

//=============================================================================
// Vector operations
//=============================================================================

Point subtract(
    const Point& a,
    const Point& b)
{
    return
    {
        a[0] - b[0],
        a[1] - b[1],
        a[2] - b[2]
    };
}


double dot(
    const Point& a,
    const Point& b)
{
    return
        a[0] * b[0] +
        a[1] * b[1] +
        a[2] * b[2];
}


Point cross(
    const Point& a,
    const Point& b)
{
    return
    {
        a[1] * b[2] -
            a[2] * b[1],

        a[2] * b[0] -
            a[0] * b[2],

        a[0] * b[1] -
            a[1] * b[0]
    };
}


double norm(
    const Point& a)
{
    return
        std::sqrt(
            dot(
                a,
                a));
}


//=============================================================================
// Triangle solid angle
//=============================================================================

double triangleSolidAngle(
    const Point& point,
    const Point& vertex0,
    const Point& vertex1,
    const Point& vertex2)
{
    const Point a =
        subtract(
            vertex0,
            point);

    const Point b =
        subtract(
            vertex1,
            point);

    const Point c =
        subtract(
            vertex2,
            point);


    const double lengthA =
        norm(
            a);

    const double lengthB =
        norm(
            b);

    const double lengthC =
        norm(
            c);


    const double numerator =
        dot(
            a,
            cross(
                b,
                c));


    const double denominator =
        lengthA *
            lengthB *
            lengthC +
        dot(
            a,
            b) *
            lengthC +
        dot(
            b,
            c) *
            lengthA +
        dot(
            c,
            a) *
            lengthB;


    return
        2.0 *
        std::atan2(
            numerator,
            denominator);
}


//=============================================================================
// Point in component
//=============================================================================

bool pointInComponent(
    const STLGeometry& geometry,
    const std::size_t componentID,
    const Point& point)
{
    const auto& topology =
        geometry.topology;

    const auto& component =
        topology.components[
            componentID];


    double solidAngle =
        0.0;


#pragma omp parallel for reduction(+:solidAngle) schedule(static)
    for(std::int64_t localFacetID = 0;
        localFacetID <
            static_cast<std::int64_t>(
                component.facets.size());
        ++localFacetID)
    {
        const std::size_t facetID =
            component.facets[
                static_cast<std::size_t>(
                    localFacetID)];


        const auto& vertexIDs =
            topology.geometry.facetVertexIDs[
                facetID];


        const Point& vertex0 =
            topology.geometry.vertices[
                vertexIDs[0]];

        const Point& vertex1 =
            topology.geometry.vertices[
                vertexIDs[1]];

        const Point& vertex2 =
            topology.geometry.vertices[
                vertexIDs[2]];


        solidAngle +=
            triangleSolidAngle(
                point,
                vertex0,
                vertex1,
                vertex2);
    }


    constexpr double pi =
        3.14159265358979323846;


    return
        std::abs(
            solidAngle) >
        2.0 * pi;
}


//=============================================================================
// Wet-side test
//=============================================================================

bool centerIsWet(
    const STLGeometry& geometry,
    const Point& point)
{
    const std::size_t componentCount =
        geometry.topology.components.size();


    for(std::size_t componentID = 0;
        componentID < componentCount;
        ++componentID)
    {
        const STLComponentFlow& componentFlow =
            geometry.flow[
                componentID];


        if(!componentFlow.active)
        {
            continue;
        }


        const bool inside =
            pointInComponent(
                geometry,
                componentID,
                point);


        if(componentFlow.fluidSide ==
           FluidSide::Inside)
        {
            if(!inside)
            {
                return false;
            }
        }
        else
        {
            if(inside)
            {
                return false;
            }
        }
    }


    return true;
}


//=============================================================================
// Separating-axis tests
//=============================================================================

bool separatedOnAxisClosed(
    const Point& vertex0,
    const Point& vertex1,
    const Point& vertex2,
    const Point& axis,
    const double halfGridSpacing)
{
    const double axisLengthSquared =
        dot(
            axis,
            axis);


    if(axisLengthSquared ==
       0.0)
    {
        return false;
    }


    const double projection0 =
        dot(
            vertex0,
            axis);

    const double projection1 =
        dot(
            vertex1,
            axis);

    const double projection2 =
        dot(
            vertex2,
            axis);


    const double triangleMin =
        std::min(
            projection0,
            std::min(
                projection1,
                projection2));

    const double triangleMax =
        std::max(
            projection0,
            std::max(
                projection1,
                projection2));


    const double boxRadius =
        halfGridSpacing *
        (
            std::abs(
                axis[0]) +
            std::abs(
                axis[1]) +
            std::abs(
                axis[2])
        );


    return
        triangleMax <
            -boxRadius ||
        triangleMin >
            boxRadius;
}


bool separatedOnAxisOpen(
    const Point& vertex0,
    const Point& vertex1,
    const Point& vertex2,
    const Point& axis,
    const double halfGridSpacing)
{
    const double axisLengthSquared =
        dot(
            axis,
            axis);


    if(axisLengthSquared ==
       0.0)
    {
        return false;
    }


    const double projection0 =
        dot(
            vertex0,
            axis);

    const double projection1 =
        dot(
            vertex1,
            axis);

    const double projection2 =
        dot(
            vertex2,
            axis);


    const double triangleMin =
        std::min(
            projection0,
            std::min(
                projection1,
                projection2));

    const double triangleMax =
        std::max(
            projection0,
            std::max(
                projection1,
                projection2));


    const double boxRadius =
        halfGridSpacing *
        (
            std::abs(
                axis[0]) +
            std::abs(
                axis[1]) +
            std::abs(
                axis[2])
        );


    return
        triangleMax <=
            -boxRadius ||
        triangleMin >=
            boxRadius;
}


//=============================================================================
// Triangle / closed-cell intersection
//=============================================================================

bool triangleIntersectsCell(
    const Point& point,
    const double gridSpacing,
    const Point& triangleVertex0,
    const Point& triangleVertex1,
    const Point& triangleVertex2)
{
    const double halfGridSpacing =
        0.5 *
        gridSpacing;


    const Point vertex0 =
        subtract(
            triangleVertex0,
            point);

    const Point vertex1 =
        subtract(
            triangleVertex1,
            point);

    const Point vertex2 =
        subtract(
            triangleVertex2,
            point);


    const double triangleMinX =
        std::min(
            vertex0[0],
            std::min(
                vertex1[0],
                vertex2[0]));

    const double triangleMaxX =
        std::max(
            vertex0[0],
            std::max(
                vertex1[0],
                vertex2[0]));


    if(triangleMaxX <
           -halfGridSpacing ||
       triangleMinX >
           halfGridSpacing)
    {
        return false;
    }


    const double triangleMinY =
        std::min(
            vertex0[1],
            std::min(
                vertex1[1],
                vertex2[1]));

    const double triangleMaxY =
        std::max(
            vertex0[1],
            std::max(
                vertex1[1],
                vertex2[1]));


    if(triangleMaxY <
           -halfGridSpacing ||
       triangleMinY >
           halfGridSpacing)
    {
        return false;
    }


    const double triangleMinZ =
        std::min(
            vertex0[2],
            std::min(
                vertex1[2],
                vertex2[2]));

    const double triangleMaxZ =
        std::max(
            vertex0[2],
            std::max(
                vertex1[2],
                vertex2[2]));


    if(triangleMaxZ <
           -halfGridSpacing ||
       triangleMinZ >
           halfGridSpacing)
    {
        return false;
    }


    const Point edge0 =
        subtract(
            vertex1,
            vertex0);

    const Point edge1 =
        subtract(
            vertex2,
            vertex1);

    const Point edge2 =
        subtract(
            vertex0,
            vertex2);


    const Point triangleNormal =
        cross(
            edge0,
            edge1);


    if(separatedOnAxisClosed(
           vertex0,
           vertex1,
           vertex2,
           triangleNormal,
           halfGridSpacing))
    {
        return false;
    }


    const std::array<Point,3> edges =
    {{
        edge0,
        edge1,
        edge2
    }};


    const std::array<Point,3> cellAxes =
    {{
        {
            1.0,
            0.0,
            0.0
        },
        {
            0.0,
            1.0,
            0.0
        },
        {
            0.0,
            0.0,
            1.0
        }
    }};


    for(const Point& edge :
        edges)
    {
        for(const Point& cellAxis :
            cellAxes)
        {
            const Point axis =
                cross(
                    edge,
                    cellAxis);


            if(separatedOnAxisClosed(
                   vertex0,
                   vertex1,
                   vertex2,
                   axis,
                   halfGridSpacing))
            {
                return false;
            }
        }
    }


    return true;
}


//=============================================================================
// Triangle / open-cell-interior intersection
//=============================================================================

bool triangleIntersectsCellInterior(
    const Point& point,
    const double gridSpacing,
    const Point& triangleVertex0,
    const Point& triangleVertex1,
    const Point& triangleVertex2)
{
    const double halfGridSpacing =
        0.5 *
        gridSpacing;


    const Point vertex0 =
        subtract(
            triangleVertex0,
            point);

    const Point vertex1 =
        subtract(
            triangleVertex1,
            point);

    const Point vertex2 =
        subtract(
            triangleVertex2,
            point);


    const double triangleMinX =
        std::min(
            vertex0[0],
            std::min(
                vertex1[0],
                vertex2[0]));

    const double triangleMaxX =
        std::max(
            vertex0[0],
            std::max(
                vertex1[0],
                vertex2[0]));


    if(triangleMaxX <=
           -halfGridSpacing ||
       triangleMinX >=
           halfGridSpacing)
    {
        return false;
    }


    const double triangleMinY =
        std::min(
            vertex0[1],
            std::min(
                vertex1[1],
                vertex2[1]));

    const double triangleMaxY =
        std::max(
            vertex0[1],
            std::max(
                vertex1[1],
                vertex2[1]));


    if(triangleMaxY <=
           -halfGridSpacing ||
       triangleMinY >=
           halfGridSpacing)
    {
        return false;
    }


    const double triangleMinZ =
        std::min(
            vertex0[2],
            std::min(
                vertex1[2],
                vertex2[2]));

    const double triangleMaxZ =
        std::max(
            vertex0[2],
            std::max(
                vertex1[2],
                vertex2[2]));


    if(triangleMaxZ <=
           -halfGridSpacing ||
       triangleMinZ >=
           halfGridSpacing)
    {
        return false;
    }


    const Point edge0 =
        subtract(
            vertex1,
            vertex0);

    const Point edge1 =
        subtract(
            vertex2,
            vertex1);

    const Point edge2 =
        subtract(
            vertex0,
            vertex2);


    const Point triangleNormal =
        cross(
            edge0,
            edge1);


    if(separatedOnAxisOpen(
           vertex0,
           vertex1,
           vertex2,
           triangleNormal,
           halfGridSpacing))
    {
        return false;
    }


    const std::array<Point,3> edges =
    {{
        edge0,
        edge1,
        edge2
    }};


    const std::array<Point,3> cellAxes =
    {{
        {
            1.0,
            0.0,
            0.0
        },
        {
            0.0,
            1.0,
            0.0
        },
        {
            0.0,
            0.0,
            1.0
        }
    }};


    for(const Point& edge :
        edges)
    {
        for(const Point& cellAxis :
            cellAxes)
        {
            const Point axis =
                cross(
                    edge,
                    cellAxis);


            if(separatedOnAxisOpen(
                   vertex0,
                   vertex1,
                   vertex2,
                   axis,
                   halfGridSpacing))
            {
                return false;
            }
        }
    }


    return true;
}


//=============================================================================
// Active surface / closed-cell intersection
//=============================================================================

bool activeSurfaceIntersectsCell(
    const STLGeometry& geometry,
    const Point& point,
    const double gridSpacing)
{
    const auto& topology =
        geometry.topology;

    const double halfGridSpacing =
        0.5 *
        gridSpacing;


    for(std::size_t componentID = 0;
        componentID <
            topology.components.size();
        ++componentID)
    {
        if(!geometry.flow[
                componentID].active)
        {
            continue;
        }


        const auto& component =
            topology.components[
                componentID];


        if(component.bounds.max[0] <
               point[0] -
                   halfGridSpacing ||
           component.bounds.min[0] >
               point[0] +
                   halfGridSpacing ||
           component.bounds.max[1] <
               point[1] -
                   halfGridSpacing ||
           component.bounds.min[1] >
               point[1] +
                   halfGridSpacing ||
           component.bounds.max[2] <
               point[2] -
                   halfGridSpacing ||
           component.bounds.min[2] >
               point[2] +
                   halfGridSpacing)
        {
            continue;
        }


        for(const std::size_t facetID :
            component.facets)
        {
            const auto& vertexIDs =
                topology.geometry.facetVertexIDs[
                    facetID];


            const Point& vertex0 =
                topology.geometry.vertices[
                    vertexIDs[0]];

            const Point& vertex1 =
                topology.geometry.vertices[
                    vertexIDs[1]];

            const Point& vertex2 =
                topology.geometry.vertices[
                    vertexIDs[2]];


            if(triangleIntersectsCell(
                   point,
                   gridSpacing,
                   vertex0,
                   vertex1,
                   vertex2))
            {
                return true;
            }
        }
    }


    return false;
}


//=============================================================================
// Active surface / open-cell-interior intersection
//=============================================================================

bool activeSurfaceIntersectsCellInterior(
    const STLGeometry& geometry,
    const Point& point,
    const double gridSpacing)
{
    const auto& topology =
        geometry.topology;

    const double halfGridSpacing =
        0.5 *
        gridSpacing;


    for(std::size_t componentID = 0;
        componentID <
            topology.components.size();
        ++componentID)
    {
        if(!geometry.flow[
                componentID].active)
        {
            continue;
        }


        const auto& component =
            topology.components[
                componentID];


        if(component.bounds.max[0] <=
               point[0] -
                   halfGridSpacing ||
           component.bounds.min[0] >=
               point[0] +
                   halfGridSpacing ||
           component.bounds.max[1] <=
               point[1] -
                   halfGridSpacing ||
           component.bounds.min[1] >=
               point[1] +
                   halfGridSpacing ||
           component.bounds.max[2] <=
               point[2] -
                   halfGridSpacing ||
           component.bounds.min[2] >=
               point[2] +
                   halfGridSpacing)
        {
            continue;
        }


        for(const std::size_t facetID :
            component.facets)
        {
            const auto& vertexIDs =
                topology.geometry.facetVertexIDs[
                    facetID];


            const Point& vertex0 =
                topology.geometry.vertices[
                    vertexIDs[0]];

            const Point& vertex1 =
                topology.geometry.vertices[
                    vertexIDs[1]];

            const Point& vertex2 =
                topology.geometry.vertices[
                    vertexIDs[2]];


            if(triangleIntersectsCellInterior(
                   point,
                   gridSpacing,
                   vertex0,
                   vertex1,
                   vertex2))
            {
                return true;
            }
        }
    }


    return false;
}

} // namespace


//=============================================================================
// STL computational-domain cell classification
//=============================================================================

CellType STLGeometry::contains(
    const Point& point,
    const double gridSpacing) const
{
    if(gridSpacing <=
       0.0)
    {
        throw std::runtime_error(
            "STLGeometry::contains requires "
            "positive grid spacing.");
    }


    if(flow.size() !=
       topology.components.size())
    {
        throw std::runtime_error(
            "STLGeometry::contains requires "
            "interpreted flow data.");
    }


    const bool wet =
        centerIsWet(
            *this,
            point);


    if(wet)
    {
        if(activeSurfaceIntersectsCell(
               *this,
               point,
               gridSpacing))
        {
            return
                CellType::Boundary;
        }


        return
            CellType::Interior;
    }


    if(activeSurfaceIntersectsCellInterior(
           *this,
           point,
           gridSpacing))
    {
        return
            CellType::Boundary;
    }


    return
        CellType::Dry;
}

} // namespace ntic::lbm::geometry