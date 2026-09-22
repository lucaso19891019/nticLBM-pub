#include "stl_geometry.hpp"
#include "stl_reader.hpp"
#include "stl_validator.hpp"

#include <cstddef>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>


namespace
{

using ntic::lbm::geometry::BoundingBox;
using ntic::lbm::geometry::Point;
using ntic::lbm::geometry::STLGeometry;


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
// Reference bounding box
//=============================================================================

BoundingBox makeReferenceBounds(
    const ntic::lbm::stl::FacetTopology& topology)
{
    BoundingBox bounds;


    for(const auto& vertex :
        topology.geometry.vertices)
    {
        bounds.expand(
            vertex);
    }


    return bounds;
}


//=============================================================================
// STL geometry construction
//=============================================================================

void testConstruction(
    const std::string& stlFile,
    const std::string& validationMode)
{
    using namespace ntic::lbm;


    //-------------------------------------------------------------------------
    // Read and validate STL
    //-------------------------------------------------------------------------

    stl::STLData data =
        stl::read(
            stlFile);


    stl::FacetTopology topology;


    stl::validate(
        data,
        validationMode,
        topology);


    require(
        !topology.geometry.vertices.empty(),
        "Validated STL topology contains no vertices.");


    //-------------------------------------------------------------------------
    // Reference data before ownership transfer
    //-------------------------------------------------------------------------

    const BoundingBox referenceBounds =
        makeReferenceBounds(
            topology);


    const std::size_t vertexCount =
        topology.geometry.vertices.size();

    const std::size_t facetCount =
        topology.geometry.facetVertexIDs.size();

    const std::size_t facetGeometryCount =
        topology.facetGeometry.size();

    const std::size_t edgeCount =
        topology.edges.size();

    const std::size_t componentIDCount =
        topology.facetComponentIDs.size();

    const std::size_t componentCount =
        topology.components.size();


    const std::size_t middleVertexIndex =
        vertexCount / 2;


    const Point firstVertex =
        topology.geometry.vertices.front();

    const Point middleVertex =
        topology.geometry.vertices[
            middleVertexIndex];

    const Point lastVertex =
        topology.geometry.vertices.back();


    //-------------------------------------------------------------------------
    // Transfer ownership to Geometry
    //-------------------------------------------------------------------------

    STLGeometry geometry =
        geometry::constructSTLGeometry(
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

    require(
        geometry.topology.facetGeometry.size() ==
            facetGeometryCount,
        "Facet geometry count changed during STL geometry construction.");

    require(
        geometry.topology.edges.size() ==
            edgeCount,
        "Edge count changed during STL geometry construction.");

    require(
        geometry.topology.facetComponentIDs.size() ==
            componentIDCount,
        "Facet component ID count changed during STL geometry construction.");

    require(
        geometry.topology.components.size() ==
            componentCount,
        "Component count changed during STL geometry construction.");


    //-------------------------------------------------------------------------
    // Coordinate preservation
    //-------------------------------------------------------------------------

    requirePointEqual(
        geometry.topology.geometry.vertices.front(),
        firstVertex,
        "First vertex changed during STL geometry construction.");

    requirePointEqual(
        geometry.topology.geometry.vertices[
            middleVertexIndex],
        middleVertex,
        "Middle vertex changed during STL geometry construction.");

    requirePointEqual(
        geometry.topology.geometry.vertices.back(),
        lastVertex,
        "Last vertex changed during STL geometry construction.");


    //-------------------------------------------------------------------------
    // Bounding box
    //-------------------------------------------------------------------------

    requirePointEqual(
        geometry.bounds.min,
        referenceBounds.min,
        "STL geometry minimum bound is incorrect.");

    requirePointEqual(
        geometry.bounds.max,
        referenceBounds.max,
        "STL geometry maximum bound is incorrect.");


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
    ntic::lbm::stl::FacetTopology topology;


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

int main(
    const int argc,
    char* argv[])
{
    if(argc != 3)
    {
        std::cerr
            << "Usage: "
            << argv[0]
            << " <stl_file> <validation_mode>"
            << std::endl;

        return 1;
    }


    try
    {
        const std::string stlFile =
            argv[1];

        const std::string validationMode =
            argv[2];


        testConstruction(
            stlFile,
            validationMode);


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