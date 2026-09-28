#include "bounding_box.hpp"
#include "flow_type.hpp"
#include "point.hpp"
#include "stl_geometry.hpp"
#include "stl_reader.hpp"
#include "stl_validator.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>


namespace
{

using ntic::lbm::geometry::BoundingBox;
using ntic::lbm::geometry::FlowType;
using ntic::lbm::geometry::Point;
using ntic::lbm::geometry::STLGeometry;


//=============================================================================
// Tolerance
//=============================================================================

constexpr double tolerance =
    1.0e-12;


//=============================================================================
// Basic checks
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


bool nearlyEqual(
    const double a,
    const double b)
{
    const double scale =
        1.0 +
        std::max(
            std::abs(a),
            std::abs(b));

    return
        std::abs(a - b) <=
        tolerance * scale;
}


bool samePoint(
    const Point& a,
    const Point& b)
{
    return
        nearlyEqual(a[0], b[0]) &&
        nearlyEqual(a[1], b[1]) &&
        nearlyEqual(a[2], b[2]);
}


Point add(
    const Point& point,
    const Point& displacement)
{
    return
    {
        point[0] + displacement[0],
        point[1] + displacement[1],
        point[2] + displacement[2]
    };
}


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


BoundingBox makeOpenBox(
    const ntic::lbm::stl::FacetTopology& topology)
{
    BoundingBox openBox;


    for(const Point& point :
        topology.geometry.vertices)
    {
        openBox.expand(
            point);
    }


    const Point padding =
    {
        5.0,
        7.0,
        9.0
    };


    openBox.min[0] -= padding[0];
    openBox.min[1] -= padding[1];
    openBox.min[2] -= padding[2];

    openBox.max[0] += padding[0];
    openBox.max[1] += padding[1];
    openBox.max[2] += padding[2];


    return openBox;
}


//=============================================================================
// Snapshot
//=============================================================================

struct TranslationSnapshot
{
    std::vector<Point>
        vertices;

    std::vector<Point>
        facetCentroids;

    std::vector<std::array<double,3>>
        facetNormals;

    std::vector<double>
        facetAreas;

    std::vector<double>
        facetQualities;

    std::vector<double>
        facetCentroidRadii;

    std::vector<std::array<std::size_t,3>>
        facetVertexIDs;

    std::vector<std::array<std::size_t,3>>
        adjacency;

    std::vector<ntic::lbm::stl::STLEdge>
        edges;

    std::vector<std::size_t>
        facetComponentIDs;

    std::vector<std::vector<std::size_t>>
        componentFacets;

    std::vector<Point>
        componentBoundsMin;

    std::vector<Point>
        componentBoundsMax;

    std::vector<Point>
        componentBoundsCenter;

    std::vector<double>
        componentBoundsScale;

    std::vector<double>
        componentSignedVolumes;

    BoundingBox bounds;

    std::vector<ntic::lbm::geometry::STLComponentContainment>
        containment;

    std::vector<ntic::lbm::geometry::STLComponentFlow>
        flow;

    FlowType flowType;
};


TranslationSnapshot captureSnapshot(
    const STLGeometry& geometry)
{
    TranslationSnapshot snapshot;


    snapshot.vertices =
        geometry.topology.geometry.vertices;


    const std::size_t facetCount =
        geometry.topology.facetGeometry.size();


    snapshot.facetCentroids.resize(
        facetCount);

    snapshot.facetNormals.resize(
        facetCount);

    snapshot.facetAreas.resize(
        facetCount);

    snapshot.facetQualities.resize(
        facetCount);

    snapshot.facetCentroidRadii.resize(
        facetCount);


    for(std::size_t facetID = 0;
        facetID < facetCount;
        ++facetID)
    {
        const auto& facet =
            geometry.topology.facetGeometry[
                facetID];


        snapshot.facetCentroids[
            facetID] =
            facet.centroid;

        snapshot.facetNormals[
            facetID] =
            facet.normal;

        snapshot.facetAreas[
            facetID] =
            facet.area;

        snapshot.facetQualities[
            facetID] =
            facet.quality;

        snapshot.facetCentroidRadii[
            facetID] =
            facet.centroidRadius;
    }


    snapshot.facetVertexIDs =
        geometry.topology.geometry.facetVertexIDs;

    snapshot.adjacency =
        geometry.topology.adjacency;

    snapshot.edges =
        geometry.topology.edges;

    snapshot.facetComponentIDs =
        geometry.topology.facetComponentIDs;


    const std::size_t componentCount =
        geometry.topology.components.size();


    snapshot.componentFacets.resize(
        componentCount);

    snapshot.componentBoundsMin.resize(
        componentCount);

    snapshot.componentBoundsMax.resize(
        componentCount);

    snapshot.componentBoundsCenter.resize(
        componentCount);

    snapshot.componentBoundsScale.resize(
        componentCount);

    snapshot.componentSignedVolumes.resize(
        componentCount);


    for(std::size_t componentID = 0;
        componentID < componentCount;
        ++componentID)
    {
        const auto& component =
            geometry.topology.components[
                componentID];


        snapshot.componentFacets[
            componentID] =
            component.facets;

        snapshot.componentBoundsMin[
            componentID] =
            component.bounds.min;

        snapshot.componentBoundsMax[
            componentID] =
            component.bounds.max;

        snapshot.componentBoundsCenter[
            componentID] =
            component.bounds.center;

        snapshot.componentBoundsScale[
            componentID] =
            component.bounds.scale;

        snapshot.componentSignedVolumes[
            componentID] =
            component.signedVolume;
    }


    snapshot.bounds =
        geometry.bounds;

    snapshot.containment =
        geometry.containment;

    snapshot.flow =
        geometry.flow;

    snapshot.flowType =
        geometry.flowType;


    return snapshot;
}


//=============================================================================
// Translation-dependent data
//=============================================================================

void checkTranslatedData(
    const STLGeometry& geometry,
    const TranslationSnapshot& before,
    const Point& displacement)
{
    require(
        geometry.topology.geometry.vertices.size() ==
            before.vertices.size(),
        "Vertex count changed during translation.");


    for(std::size_t vertexID = 0;
        vertexID < before.vertices.size();
        ++vertexID)
    {
        require(
            samePoint(
                geometry.topology.geometry.vertices[
                    vertexID],
                add(
                    before.vertices[
                        vertexID],
                    displacement)),
            "Canonical vertex was not translated correctly.");
    }


    require(
        geometry.topology.facetGeometry.size() ==
            before.facetCentroids.size(),
        "Facet count changed during translation.");


    for(std::size_t facetID = 0;
        facetID < before.facetCentroids.size();
        ++facetID)
    {
        require(
            samePoint(
                geometry.topology.facetGeometry[
                    facetID].centroid,
                add(
                    before.facetCentroids[
                        facetID],
                    displacement)),
            "Facet centroid was not translated correctly.");
    }


    require(
        geometry.topology.components.size() ==
            before.componentBoundsMin.size(),
        "Component count changed during translation.");


    for(std::size_t componentID = 0;
        componentID <
            before.componentBoundsMin.size();
        ++componentID)
    {
        const auto& bounds =
            geometry.topology.components[
                componentID].bounds;


        require(
            samePoint(
                bounds.min,
                add(
                    before.componentBoundsMin[
                        componentID],
                    displacement)),
            "Component bounds.min was not translated correctly.");


        require(
            samePoint(
                bounds.max,
                add(
                    before.componentBoundsMax[
                        componentID],
                    displacement)),
            "Component bounds.max was not translated correctly.");


        require(
            samePoint(
                bounds.center,
                add(
                    before.componentBoundsCenter[
                        componentID],
                    displacement)),
            "Component bounds.center was not translated correctly.");
    }


    require(
        samePoint(
            geometry.bounds.min,
            add(
                before.bounds.min,
                displacement)),
        "STL geometry bounds.min was not translated correctly.");


    require(
        samePoint(
            geometry.bounds.max,
            add(
                before.bounds.max,
                displacement)),
        "STL geometry bounds.max was not translated correctly.");
}


//=============================================================================
// Translation-invariant data
//=============================================================================

void checkInvariantData(
    const STLGeometry& geometry,
    const TranslationSnapshot& before)
{
    require(
        geometry.topology.geometry.facetVertexIDs ==
            before.facetVertexIDs,
        "Facet vertex IDs changed during translation.");


    require(
        geometry.topology.adjacency ==
            before.adjacency,
        "Facet adjacency changed during translation.");


    require(
        geometry.topology.facetComponentIDs ==
            before.facetComponentIDs,
        "Facet component IDs changed during translation.");


    require(
        geometry.topology.edges.size() ==
            before.edges.size(),
        "STL edge count changed during translation.");


    for(std::size_t edgeID = 0;
        edgeID < before.edges.size();
        ++edgeID)
    {
        const auto& oldEdge =
            before.edges[
                edgeID];

        const auto& newEdge =
            geometry.topology.edges[
                edgeID];


        require(
            newEdge.v0 ==
                oldEdge.v0 &&
            newEdge.v1 ==
                oldEdge.v1 &&
            newEdge.facets ==
                oldEdge.facets,
            "STL edge topology changed during translation.");
    }


    for(std::size_t facetID = 0;
        facetID < before.facetNormals.size();
        ++facetID)
    {
        const auto& facet =
            geometry.topology.facetGeometry[
                facetID];


        require(
            facet.normal ==
                before.facetNormals[
                    facetID],
            "Facet normal changed during translation.");


        require(
            facet.area ==
                before.facetAreas[
                    facetID],
            "Facet area changed during translation.");


        require(
            facet.quality ==
                before.facetQualities[
                    facetID],
            "Facet quality changed during translation.");


        require(
            facet.centroidRadius ==
                before.facetCentroidRadii[
                    facetID],
            "Facet centroid radius changed during translation.");
    }


    for(std::size_t componentID = 0;
        componentID <
            before.componentFacets.size();
        ++componentID)
    {
        const auto& component =
            geometry.topology.components[
                componentID];


        require(
            component.facets ==
                before.componentFacets[
                    componentID],
            "Component facet IDs changed during translation.");


        require(
            component.bounds.scale ==
                before.componentBoundsScale[
                    componentID],
            "Component bounds scale changed during translation.");


        require(
            component.signedVolume ==
                before.componentSignedVolumes[
                    componentID],
            "Component signed volume changed during translation.");
    }


    require(
        geometry.containment.size() ==
            before.containment.size(),
        "Containment size changed during translation.");


    for(std::size_t componentID = 0;
        componentID <
            before.containment.size();
        ++componentID)
    {
        const auto& oldContainment =
            before.containment[
                componentID];

        const auto& newContainment =
            geometry.containment[
                componentID];


        require(
            newContainment.parent ==
                oldContainment.parent &&
            newContainment.children ==
                oldContainment.children &&
            newContainment.root ==
                oldContainment.root &&
            newContainment.level ==
                oldContainment.level,
            "Containment data changed during translation.");
    }


    require(
        geometry.flow.size() ==
            before.flow.size(),
        "Flow data size changed during translation.");


    for(std::size_t componentID = 0;
        componentID <
            before.flow.size();
        ++componentID)
    {
        require(
            geometry.flow[
                componentID].active ==
                before.flow[
                    componentID].active &&
            geometry.flow[
                componentID].fluidSide ==
                before.flow[
                    componentID].fluidSide,
            "Flow data changed during translation.");
    }


    require(
        geometry.flowType ==
            before.flowType,
        "Flow type changed during translation.");
}


//=============================================================================
// Internal translation
//=============================================================================

void testInternal(
    STLGeometry& geometry)
{
    require(
        geometry.flowType ==
            FlowType::Internal,
        "Internal translation requires internal STL geometry.");


    const Point targetPoint =
    {
        2.5,
        4.0,
        6.5
    };


    const TranslationSnapshot before =
        captureSnapshot(
            geometry);


    const Point displacement =
        subtract(
            targetPoint,
            before.bounds.min);


    geometry.translate(
        &targetPoint);


    checkTranslatedData(
        geometry,
        before,
        displacement);


    checkInvariantData(
        geometry,
        before);


    require(
        samePoint(
            geometry.bounds.min,
            targetPoint),
        "Internal translation did not place "
        "geometry.bounds.min at targetPoint.");
}


//=============================================================================
// External translation
//=============================================================================

void testExternal(
    STLGeometry& geometry)
{
    require(
        geometry.flowType ==
            FlowType::External,
        "External translation requires external STL geometry.");


    const TranslationSnapshot before =
        captureSnapshot(
            geometry);


    const BoundingBox openBoxBefore =
        geometry.openBox;


    const Point relativeMinBefore =
        subtract(
            geometry.bounds.min,
            geometry.openBox.min);


    const Point relativeMaxBefore =
        subtract(
            geometry.bounds.max,
            geometry.openBox.min);


    const Point displacement =
    {
        -geometry.openBox.min[0],
        -geometry.openBox.min[1],
        -geometry.openBox.min[2]
    };


    geometry.translate();


    checkTranslatedData(
        geometry,
        before,
        displacement);


    checkInvariantData(
        geometry,
        before);


    const Point origin =
    {
        0.0,
        0.0,
        0.0
    };


    require(
        samePoint(
            geometry.openBox.min,
            origin),
        "External translation did not place "
        "geometry.openBox.min at the origin.");


    require(
        samePoint(
            geometry.openBox.max,
            add(
                openBoxBefore.max,
                displacement)),
        "External openBox.max was not translated correctly.");


    const Point relativeMinAfter =
        subtract(
            geometry.bounds.min,
            geometry.openBox.min);


    const Point relativeMaxAfter =
        subtract(
            geometry.bounds.max,
            geometry.openBox.min);


    require(
        samePoint(
            relativeMinAfter,
            relativeMinBefore),
        "STL relative position to openBox changed "
        "during external translation.");


    require(
        samePoint(
            relativeMaxAfter,
            relativeMaxBefore),
        "STL relative extent to openBox changed "
        "during external translation.");
}


//=============================================================================
// Repeated translation
//=============================================================================

void testRepeated(
    STLGeometry& geometry)
{
    require(
        geometry.flowType ==
            FlowType::Internal,
        "Repeated translation test requires internal STL geometry.");


    const Point targetPoint0 =
    {
        1.25,
        -2.50,
        3.75
    };


    const Point targetPoint1 =
    {
        -4.50,
        5.25,
        -6.00
    };


    const TranslationSnapshot before =
        captureSnapshot(
            geometry);


    const Point totalDisplacement =
        subtract(
            targetPoint1,
            before.bounds.min);


    geometry.translate(
        &targetPoint0);


    geometry.translate(
        &targetPoint1);


    checkTranslatedData(
        geometry,
        before,
        totalDisplacement);


    checkInvariantData(
        geometry,
        before);


    require(
        samePoint(
            geometry.bounds.min,
            targetPoint1),
        "Repeated internal translation did not place "
        "geometry.bounds.min at the final target point.");
}


//=============================================================================
// Test dispatch
//=============================================================================

void runTest(
    STLGeometry& geometry,
    const std::string& testCase)
{


    if(testCase ==
       "internal")
    {
        testInternal(
            geometry);

        return;
    }


    if(testCase ==
       "external")
    {
        testExternal(
            geometry);

        return;
    }


    if(testCase ==
       "repeated")
    {
        testRepeated(
            geometry);

        return;
    }


    throw std::runtime_error(
        "Unknown STL translation test case: " +
        testCase);
}

} // namespace


//=============================================================================
// Main
//=============================================================================

int main(
    const int argc,
    char** argv)
{
    if(argc != 3)
    {
        std::cerr
            << "Usage: "
            << argv[0]
            << " <test_case> <stl_file>"
            << std::endl;

        return 1;
    }


    const std::string testCase =
        argv[1];

    const std::string stlFile =
        argv[2];


    try
    {
        ntic::lbm::stl::STLData data =
            ntic::lbm::stl::read(
                stlFile);


        ntic::lbm::stl::FacetTopology topology;


        ntic::lbm::stl::validate(
            data,
            "full",
            topology);


        FlowType flowType =
            FlowType::Internal;


        if(testCase ==
           "external")
        {
            flowType =
                FlowType::External;
        }


        BoundingBox openBox;

        const BoundingBox* openBoxPointer =
            nullptr;


        if(flowType ==
           FlowType::External)
        {
            openBox =
                makeOpenBox(
                    topology);

            openBoxPointer =
                &openBox;
        }


        STLGeometry geometry(
            std::move(topology),
            flowType,
            openBoxPointer);


        runTest(
            geometry,
            testCase);


        std::cout
            << "[PASS] "
            << testCase
            << " : "
            << stlFile
            << std::endl;
    }
    catch(const std::exception& error)
    {
        std::cerr
            << "[FAIL] "
            << testCase
            << " : "
            << stlFile
            << std::endl
            << "       "
            << error.what()
            << std::endl;

        return 1;
    }


    return 0;
}