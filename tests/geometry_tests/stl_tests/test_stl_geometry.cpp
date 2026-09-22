#include "stl_geometry.hpp"

#include <cstddef>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>


namespace
{

using ntic::lbm::geometry::Point;
using ntic::lbm::geometry::STLGeometry;
using ntic::lbm::stl::FacetTopology;


//=============================================================================
// Test utilities
//=============================================================================

void require(
    const bool condition,
    const std::string& message)
{
    if(!condition)
    {
        throw std::runtime_error(
            message);
    }
}


void requirePointEqual(
    const Point& actual,
    const Point& expected,
    const std::string& message)
{
    for(std::size_t d = 0;
        d < 3;
        ++d)
    {
        if(actual[d] !=
           expected[d])
        {
            throw std::runtime_error(
                message);
        }
    }
}


//=============================================================================
// Test topology
//=============================================================================

FacetTopology makeTestTopology()
{
    FacetTopology topology;


    topology.geometry.vertices = {

        {-2.0, 3.0, 5.0},
        { 4.0, 3.0, 5.0},
        {-2.0, 8.0, 5.0},
        { 4.0, 8.0, 5.0},

        {-2.0, 3.0, 12.0},
        { 4.0, 3.0, 12.0},
        {-2.0, 8.0, 12.0},
        { 4.0, 8.0, 12.0}
    };


    topology.geometry.facetVertexIDs = {

        {0, 3, 1},
        {0, 2, 3},

        {4, 5, 7},
        {4, 7, 6},

        {0, 1, 5},
        {0, 5, 4},

        {2, 6, 7},
        {2, 7, 3},

        {0, 4, 6},
        {0, 6, 2},

        {1, 3, 7},
        {1, 7, 5}
    };


    return topology;
}


//=============================================================================
// STL geometry construction
//=============================================================================

void testConstruction()
{
    FacetTopology topology =
        makeTestTopology();


    const std::size_t vertexCount =
        topology.geometry.vertices.size();

    const std::size_t facetCount =
        topology.geometry.facetVertexIDs.size();


    const Point firstVertex =
        topology.geometry.vertices.front();

    const std::size_t middleIndex =
        vertexCount / 2;

    const Point middleVertex =
        topology.geometry.vertices[
            middleIndex];

    const Point lastVertex =
        topology.geometry.vertices.back();


    STLGeometry geometry =
        ntic::lbm::geometry::constructSTLGeometry(
            std::move(topology));


    //-------------------------------------------------------------------------
    // Topology preservation
    //-------------------------------------------------------------------------

    require(
        geometry.topology.geometry.vertices.size() ==
            vertexCount,
        "Vertex count changed during STL geometry construction.");

    require(
        geometry.topology.geometry.facetVertexIDs.size() ==
            facetCount,
        "Facet count changed during STL geometry construction.");


    //-------------------------------------------------------------------------
    // Coordinate preservation
    //-------------------------------------------------------------------------

    requirePointEqual(
        geometry.topology.geometry.vertices.front(),
        firstVertex,
        "First vertex changed during STL geometry construction.");

    requirePointEqual(
        geometry.topology.geometry.vertices[
            middleIndex],
        middleVertex,
        "Middle vertex changed during STL geometry construction.");

    requirePointEqual(
        geometry.topology.geometry.vertices.back(),
        lastVertex,
        "Last vertex changed during STL geometry construction.");


    //-------------------------------------------------------------------------
    // Tight bounding box
    //-------------------------------------------------------------------------

    const Point expectedMin{
        -2.0,
        3.0,
        5.0
    };

    const Point expectedMax{
        4.0,
        8.0,
        12.0
    };


    requirePointEqual(
        geometry.bounds.min,
        expectedMin,
        "STL geometry minimum bound is incorrect.");

    requirePointEqual(
        geometry.bounds.max,
        expectedMax,
        "STL geometry maximum bound is incorrect.");


    //-------------------------------------------------------------------------
    // Bounding-box containment
    //-------------------------------------------------------------------------

    for(const auto& vertex :
        geometry.topology.geometry.vertices)
    {
        require(
            geometry.bounds.contains(
                vertex),
            "STL geometry bounds do not contain all vertices.");
    }
}


//=============================================================================
// Empty STL geometry
//=============================================================================

void testEmptyGeometry()
{
    FacetTopology topology;


    bool caughtExpectedException =
        false;


    try
    {
        STLGeometry geometry =
            ntic::lbm::geometry::constructSTLGeometry(
                std::move(topology));

        static_cast<void>(
            geometry);
    }
    catch(const std::runtime_error&)
    {
        caughtExpectedException =
            true;
    }


    require(
        caughtExpectedException,
        "Empty STL geometry construction did not throw.");
}

} // namespace


//=============================================================================
// Main
//=============================================================================

int main()
{
    try
    {
        testConstruction();

        testEmptyGeometry();


        std::cout
            << "STL geometry tests passed."
            << std::endl;


        return 0;
    }
    catch(const std::exception& exception)
    {
        std::cerr
            << "STL geometry tests failed: "
            << exception.what()
            << std::endl;


        return 1;
    }
}