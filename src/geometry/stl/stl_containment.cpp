#include "stl_containment.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>


namespace ntic::lbm::geometry
{

namespace
{

//=============================================================================
// Vector operations
//=============================================================================

std::array<double,3> subtract(
    const std::array<double,3>& a,
    const std::array<double,3>& b)
{
    return
    {{
        a[0] - b[0],
        a[1] - b[1],
        a[2] - b[2]
    }};
}


double dot(
    const std::array<double,3>& a,
    const std::array<double,3>& b)
{
    return
        a[0] * b[0] +
        a[1] * b[1] +
        a[2] * b[2];
}


std::array<double,3> cross(
    const std::array<double,3>& a,
    const std::array<double,3>& b)
{
    return
    {{
        a[1] * b[2] -
            a[2] * b[1],

        a[2] * b[0] -
            a[0] * b[2],

        a[0] * b[1] -
            a[1] * b[0]
    }};
}


double norm(
    const std::array<double,3>& a)
{
    return std::sqrt(
        dot(
            a,
            a));
}


//=============================================================================
// Component bounding-box containment
//=============================================================================

bool boundsContain(
    const stl::GeometryBounds& outer,
    const stl::GeometryBounds& inner)
{
    return
        inner.min[0] >= outer.min[0] &&
        inner.min[1] >= outer.min[1] &&
        inner.min[2] >= outer.min[2] &&

        inner.max[0] <= outer.max[0] &&
        inner.max[1] <= outer.max[1] &&
        inner.max[2] <= outer.max[2];
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
    const auto a =
        subtract(
            vertex0,
            point);

    const auto b =
        subtract(
            vertex1,
            point);

    const auto c =
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
// Point in closed component
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


        solidAngle +=
            triangleSolidAngle(
                point,
                vertex0,
                vertex1,
                vertex2);
    }


    constexpr double PI =
        3.141592653589793238462643383279502884;


    return
        std::abs(
            solidAngle) >
        2.0 * PI;
}


//=============================================================================
// Component vertex IDs
//=============================================================================

std::vector<std::size_t> collectComponentVertices(
    const STLGeometry& geometry,
    const std::size_t componentID)
{
    const auto& topology =
        geometry.topology;

    const auto& component =
        topology.components[
            componentID];


    const std::size_t vertexCount =
        topology.geometry.vertices.size();


    std::vector<unsigned char> used(
        vertexCount,
        0);


    for(const std::size_t facetID :
        component.facets)
    {
        const auto& vertexIDs =
            topology.geometry.facetVertexIDs[
                facetID];


        used[vertexIDs[0]] =
            1;

        used[vertexIDs[1]] =
            1;

        used[vertexIDs[2]] =
            1;
    }


    std::vector<std::size_t> vertexIDs;


    for(std::size_t vertexID = 0;
        vertexID < vertexCount;
        ++vertexID)
    {
        if(used[vertexID] != 0)
        {
            vertexIDs.push_back(
                vertexID);
        }
    }


    return vertexIDs;
}


//=============================================================================
// Component in component
//=============================================================================

bool componentInComponent(
    const STLGeometry& geometry,
    const std::size_t innerComponentID,
    const std::size_t outerComponentID)
{
    const std::vector<std::size_t> vertexIDs =
        collectComponentVertices(
            geometry,
            innerComponentID);


    if(vertexIDs.empty())
    {
        throw std::runtime_error(
            "STL component contains no vertices.");
    }


    std::vector<unsigned char> inside(
        vertexIDs.size(),
        0);


#pragma omp parallel for schedule(static)
    for(std::ptrdiff_t index = 0;
        index <
            static_cast<std::ptrdiff_t>(
                vertexIDs.size());
        ++index)
    {
        const std::size_t vertexID =
            vertexIDs[
                static_cast<std::size_t>(
                    index)];


        const Point& point =
            geometry.topology.geometry.vertices[
                vertexID];


        inside[
            static_cast<std::size_t>(
                index)] =
            pointInComponent(
                geometry,
                outerComponentID,
                point)
                ? 1
                : 0;
    }


    for(const unsigned char result :
        inside)
    {
        if(result == 0)
        {
            return false;
        }
    }


    return true;
}


//=============================================================================
// Direct parent
//=============================================================================

std::size_t findDirectParent(
    const std::vector<std::vector<unsigned char>>& contains,
    const std::size_t componentID)
{
    const std::size_t componentCount =
        contains.size();


    std::size_t parent =
        std::numeric_limits<std::size_t>::max();


    for(std::size_t candidate = 0;
        candidate < componentCount;
        ++candidate)
    {
        if(candidate ==
           componentID)
        {
            continue;
        }


        if(contains[candidate][componentID] ==
           0)
        {
            continue;
        }


        bool direct =
            true;


        for(std::size_t intermediate = 0;
            intermediate < componentCount;
            ++intermediate)
        {
            if(intermediate ==
                   candidate ||
               intermediate ==
                   componentID)
            {
                continue;
            }


            if(contains[candidate][intermediate] != 0 &&
               contains[intermediate][componentID] != 0)
            {
                direct =
                    false;

                break;
            }
        }


        if(!direct)
        {
            continue;
        }


        if(parent !=
           std::numeric_limits<std::size_t>::max())
        {
            throw std::runtime_error(
                "Invalid STL component containment: "
                "component has multiple direct parents.");
        }


        parent =
            candidate;
    }


    return parent;
}


//=============================================================================
// Root and level
//=============================================================================

void assignRootAndLevel(
    std::vector<STLComponentContainment>& containment,
    const std::size_t componentID)
{
    const std::size_t noParent =
        std::numeric_limits<std::size_t>::max();


    std::size_t current =
        componentID;

    std::size_t level =
        0;


    for(std::size_t count = 0;
        count <= containment.size();
        ++count)
    {
        const std::size_t parent =
            containment[current].parent;


        if(parent ==
           noParent)
        {
            containment[componentID].root =
                current;

            containment[componentID].level =
                level;

            return;
        }


        current =
            parent;

        ++level;
    }


    throw std::runtime_error(
        "Invalid STL component containment: "
        "containment cycle detected.");
}

} // namespace


//=============================================================================
// STL component containment analysis
//=============================================================================

void analyzeSTLContainment(
    STLGeometry& geometry)
{
    const std::size_t componentCount =
        geometry.topology.components.size();


    geometry.containment.clear();

    geometry.containment.resize(
        componentCount);


    if(componentCount == 0)
    {
        throw std::runtime_error(
            "Cannot analyze STL containment without components.");
    }


    //-------------------------------------------------------------------------
    // Complete containment relation
    //-------------------------------------------------------------------------

    std::vector<std::vector<unsigned char>> contains(
        componentCount,
        std::vector<unsigned char>(
            componentCount,
            0));


    for(std::size_t outerComponentID = 0;
        outerComponentID < componentCount;
        ++outerComponentID)
    {
        const auto& outerBounds =
            geometry.topology.components[
                outerComponentID].bounds;


        for(std::size_t innerComponentID = 0;
            innerComponentID < componentCount;
            ++innerComponentID)
        {
            if(innerComponentID ==
               outerComponentID)
            {
                continue;
            }


            const auto& innerBounds =
                geometry.topology.components[
                    innerComponentID].bounds;


            if(!boundsContain(
                   outerBounds,
                   innerBounds))
            {
                continue;
            }


            if(componentInComponent(
                   geometry,
                   innerComponentID,
                   outerComponentID))
            {
                contains[
                    outerComponentID][
                    innerComponentID] =
                    1;
            }
        }
    }


    //-------------------------------------------------------------------------
    // Direct parents
    //-------------------------------------------------------------------------

    for(std::size_t componentID = 0;
        componentID < componentCount;
        ++componentID)
    {
        geometry.containment[
            componentID].parent =
            findDirectParent(
                contains,
                componentID);
    }


    //-------------------------------------------------------------------------
    // Children
    //-------------------------------------------------------------------------

    const std::size_t noParent =
        std::numeric_limits<std::size_t>::max();


    for(std::size_t componentID = 0;
        componentID < componentCount;
        ++componentID)
    {
        const std::size_t parent =
            geometry.containment[
                componentID].parent;


        if(parent !=
           noParent)
        {
            geometry.containment[
                parent].children.push_back(
                    componentID);
        }
    }


    //-------------------------------------------------------------------------
    // Roots and levels
    //-------------------------------------------------------------------------

    for(std::size_t componentID = 0;
        componentID < componentCount;
        ++componentID)
    {
        assignRootAndLevel(
            geometry.containment,
            componentID);
    }


    //-------------------------------------------------------------------------
    // Maximum supported containment depth
    //-------------------------------------------------------------------------

    for(std::size_t componentID = 0;
        componentID < componentCount;
        ++componentID)
    {
        if(geometry.containment[
               componentID].level >= 2)
        {
            throw std::runtime_error(
                "Invalid STL component containment: "
                "three or more nesting levels are not supported.");
        }
    }
}

} // namespace ntic::lbm::geometry