#include "flow_type.hpp"
#include "stl_geometry.hpp"
#include "stl_reader.hpp"
#include "stl_validator.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>


namespace
{

using ntic::lbm::geometry::BoundingBox;
using ntic::lbm::geometry::Point;

using ntic::lbm::geometry::FlowType;
using ntic::lbm::geometry::FluidSide;
using ntic::lbm::geometry::STLGeometry;




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


    const double padding =
        1.0;

    openBox.min[0] -= padding;
    openBox.min[1] -= padding;
    openBox.min[2] -= padding;

    openBox.max[0] += padding;
    openBox.max[1] += padding;
    openBox.max[2] += padding;


    return openBox;
}

//=============================================================================
// Test helpers
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


bool sameVector(
    const std::array<double,3>& a,
    const std::array<double,3>& b)
{
    return
        a[0] == b[0] &&
        a[1] == b[1] &&
        a[2] == b[2];
}


bool oppositeVector(
    const std::array<double,3>& a,
    const std::array<double,3>& b)
{
    return
        a[0] == -b[0] &&
        a[1] == -b[1] &&
        a[2] == -b[2];
}


//=============================================================================
// Root count
//=============================================================================

std::size_t countRoots(
    const STLGeometry& geometry)
{
    const std::size_t noParent =
        std::numeric_limits<std::size_t>::max();


    std::size_t rootCount =
        0;


    for(const auto& component :
        geometry.containment)
    {
        if(component.parent ==
           noParent)
        {
            ++rootCount;
        }
    }


    return rootCount;
}


//=============================================================================
// Orientation snapshot
//=============================================================================

struct OrientationSnapshot
{
    std::vector<std::array<std::size_t,3>>
        facetVertexIDs;

    std::vector<std::array<double,3>>
        normals;

    std::vector<double>
        signedVolumes;
};


OrientationSnapshot captureOrientation(
    const STLGeometry& geometry)
{
    OrientationSnapshot snapshot;


    snapshot.facetVertexIDs =
        geometry.topology.geometry.facetVertexIDs;


    snapshot.normals.resize(
        geometry.topology.facetGeometry.size());


    for(std::size_t facetID = 0;
        facetID <
            geometry.topology.facetGeometry.size();
        ++facetID)
    {
        snapshot.normals[facetID] =
            geometry.topology.facetGeometry[
                facetID].normal;
    }


    snapshot.signedVolumes.resize(
        geometry.topology.components.size());


    for(std::size_t componentID = 0;
        componentID <
            geometry.topology.components.size();
        ++componentID)
    {
        snapshot.signedVolumes[
            componentID] =
            geometry.topology.components[
                componentID].signedVolume;
    }


    return snapshot;
}


//=============================================================================
// Topology snapshot
//=============================================================================

struct TopologySnapshot
{
    std::vector<std::array<double,3>>
        vertices;

    std::vector<std::array<std::size_t,3>>
        adjacency;

    std::vector<ntic::lbm::stl::STLEdge>
        edges;

    std::vector<std::size_t>
        facetComponentIDs;
};


TopologySnapshot captureTopology(
    const STLGeometry& geometry)
{
    TopologySnapshot snapshot;


    snapshot.vertices =
        geometry.topology.geometry.vertices;

    snapshot.adjacency =
        geometry.topology.adjacency;

    snapshot.edges =
        geometry.topology.edges;

    snapshot.facetComponentIDs =
        geometry.topology.facetComponentIDs;


    return snapshot;
}


//=============================================================================
// Orientation-independent topology check
//=============================================================================

void checkTopologyUnchanged(
    const STLGeometry& geometry,
    const TopologySnapshot& snapshot)
{
    require(
        geometry.topology.geometry.vertices ==
            snapshot.vertices,
        "Canonical vertices changed during flow interpretation.");


    require(
        geometry.topology.adjacency ==
            snapshot.adjacency,
        "Facet adjacency changed during flow interpretation.");


    require(
        geometry.topology.facetComponentIDs ==
            snapshot.facetComponentIDs,
        "Facet component IDs changed during flow interpretation.");


    require(
        geometry.topology.edges.size() ==
            snapshot.edges.size(),
        "STL edge count changed during flow interpretation.");


    for(std::size_t edgeID = 0;
        edgeID < snapshot.edges.size();
        ++edgeID)
    {
        const auto& before =
            snapshot.edges[
                edgeID];

        const auto& after =
            geometry.topology.edges[
                edgeID];


        require(
            after.v0 == before.v0 &&
            after.v1 == before.v1 &&
            after.facets == before.facets,
            "STL edge topology changed during flow interpretation.");
    }
}


//=============================================================================
// Component orientation transformation check
//=============================================================================

void checkComponentOrientationChange(
    const STLGeometry& geometry,
    const OrientationSnapshot& before,
    const std::size_t componentID,
    const bool shouldFlip)
{
    const auto& component =
        geometry.topology.components[
            componentID];


    const double oldVolume =
        before.signedVolumes[
            componentID];

    const double newVolume =
        component.signedVolume;


    if(shouldFlip)
    {
        require(
            newVolume == -oldVolume,
            "Component signed volume was not reversed.");
    }
    else
    {
        require(
            newVolume == oldVolume,
            "Component signed volume changed unexpectedly.");
    }


    for(const std::size_t facetID :
        component.facets)
    {
        const auto& oldVertexIDs =
            before.facetVertexIDs[
                facetID];

        const auto& newVertexIDs =
            geometry.topology.geometry.facetVertexIDs[
                facetID];


        const auto& oldNormal =
            before.normals[
                facetID];

        const auto& newNormal =
            geometry.topology.facetGeometry[
                facetID].normal;


        if(shouldFlip)
        {
            require(
                newVertexIDs[0] ==
                    oldVertexIDs[0] &&
                newVertexIDs[1] ==
                    oldVertexIDs[2] &&
                newVertexIDs[2] ==
                    oldVertexIDs[1],
                "Facet winding was not reversed correctly.");


            require(
                oppositeVector(
                    newNormal,
                    oldNormal),
                "Facet normal was not reversed correctly.");
        }
        else
        {
            require(
                newVertexIDs ==
                    oldVertexIDs,
                "Facet winding changed unexpectedly.");


            require(
                sameVector(
                    newNormal,
                    oldNormal),
                "Facet normal changed unexpectedly.");
        }
    }
}


//=============================================================================
// Expected orientation
//=============================================================================

bool expectedOutward(
    const STLGeometry& geometry,
    const std::size_t componentID)
{
    const auto& flow =
        geometry.flow[
            componentID];


    require(
        flow.active,
        "Cannot request orientation for inactive component.");


    return
        flow.fluidSide ==
            FluidSide::Outside;
}


void checkActiveOrientations(
    const STLGeometry& geometry)
{
    for(std::size_t componentID = 0;
        componentID <
            geometry.topology.components.size();
        ++componentID)
    {
        const auto& flow =
            geometry.flow[
                componentID];


        if(!flow.active)
        {
            continue;
        }


        const double signedVolume =
            geometry.topology.components[
                componentID].signedVolume;


        require(
            signedVolume != 0.0,
            "Active component has zero signed volume.");


        const bool outward =
            signedVolume > 0.0;


        require(
            outward ==
                expectedOutward(
                    geometry,
                    componentID),
            "Active component orientation does not match "
            "its fluid side.");
    }
}


//=============================================================================
// Internal flow semantics
//=============================================================================

void checkInternalSemantics(
    const STLGeometry& geometry)
{
    require(
        countRoots(
            geometry) == 1,
        "Internal test geometry must contain exactly one root.");


    require(
        geometry.flow.size() ==
            geometry.topology.components.size(),
        "Internal flow data size is incorrect.");


    for(std::size_t componentID = 0;
        componentID <
            geometry.containment.size();
        ++componentID)
    {
        const std::size_t level =
            geometry.containment[
                componentID].level;

        const auto& flow =
            geometry.flow[
                componentID];


        require(
            flow.active,
            "Internal flow component is unexpectedly inactive.");


        if(level == 0)
        {
            require(
                flow.fluidSide ==
                    FluidSide::Inside,
                "Internal root must have fluid on the inside.");
        }
        else if(level == 1)
        {
            require(
                flow.fluidSide ==
                    FluidSide::Outside,
                "Internal child must have fluid on the outside.");
        }
        else
        {
            throw std::runtime_error(
                "Unexpected containment level in internal test.");
        }
    }


    checkActiveOrientations(
        geometry);
}


//=============================================================================
// External flow semantics
//=============================================================================

void checkExternalSemantics(
    const STLGeometry& geometry)
{
    require(
        geometry.flow.size() ==
            geometry.topology.components.size(),
        "External flow data size is incorrect.");


    for(std::size_t componentID = 0;
        componentID <
            geometry.containment.size();
        ++componentID)
    {
        const std::size_t level =
            geometry.containment[
                componentID].level;

        const auto& flow =
            geometry.flow[
                componentID];


        if(level == 0)
        {
            require(
                flow.active,
                "External root must be active.");


            require(
                flow.fluidSide ==
                    FluidSide::Outside,
                "External root must have fluid on the outside.");
        }
        else if(level == 1)
        {
            require(
                !flow.active,
                "External child must be inactive.");
        }
        else
        {
            throw std::runtime_error(
                "Unexpected containment level in external test.");
        }
    }


    checkActiveOrientations(
        geometry);
}


//=============================================================================
// Internal root
//=============================================================================

void testInternalRoot(
    STLGeometry& geometry)
{
    require(
        geometry.flowType ==
            FlowType::Internal,
        "internal_root requires internal flow type.");


    require(
        geometry.topology.components.size() == 1,
        "internal_root requires exactly one component.");


    checkInternalSemantics(
        geometry);


    const OrientationSnapshot before =
        captureOrientation(
            geometry);

    const TopologySnapshot topologyBefore =
        captureTopology(
            geometry);


    geometry.interpretFlow();


    checkInternalSemantics(
        geometry);


    checkTopologyUnchanged(
        geometry,
        topologyBefore);


    checkComponentOrientationChange(
        geometry,
        before,
        0,
        false);
}

//=============================================================================
// Internal nested
//=============================================================================

void testInternalNested(
    STLGeometry& geometry)
{
    require(
        geometry.flowType ==
            FlowType::Internal,
        "internal_nested requires internal flow type.");


    require(
        countRoots(
            geometry) == 1,
        "internal_nested requires exactly one root.");


    checkInternalSemantics(
        geometry);


    const OrientationSnapshot before =
        captureOrientation(
            geometry);

    const TopologySnapshot topologyBefore =
        captureTopology(
            geometry);


    geometry.interpretFlow();


    checkInternalSemantics(
        geometry);


    checkTopologyUnchanged(
        geometry,
        topologyBefore);


    for(std::size_t componentID = 0;
        componentID <
            geometry.topology.components.size();
        ++componentID)
    {
        checkComponentOrientationChange(
            geometry,
            before,
            componentID,
            false);
    }
}

//=============================================================================
// Internal multiple roots
//=============================================================================

void testInternalMultipleRoots(
    ntic::lbm::stl::FacetTopology topology)
{
    bool threw =
        false;


    try
    {
        STLGeometry geometry(
            std::move(topology),
            FlowType::Internal);

        static_cast<void>(
            geometry);
    }
    catch(const std::runtime_error&)
    {
        threw =
            true;
    }


    require(
        threw,
        "Internal flow accepted multiple root components.");
}

//=============================================================================
// External roots
//=============================================================================

void testExternalRoots(
    STLGeometry& geometry)
{
    require(
        geometry.flowType ==
            FlowType::External,
        "external_roots requires external flow type.");


    require(
        countRoots(
            geometry) >= 1,
        "external_roots requires at least one root.");


    checkExternalSemantics(
        geometry);


    const OrientationSnapshot before =
        captureOrientation(
            geometry);

    const TopologySnapshot topologyBefore =
        captureTopology(
            geometry);


    geometry.interpretFlow();


    checkExternalSemantics(
        geometry);


    checkTopologyUnchanged(
        geometry,
        topologyBefore);


    for(std::size_t componentID = 0;
        componentID <
            geometry.topology.components.size();
        ++componentID)
    {
        checkComponentOrientationChange(
            geometry,
            before,
            componentID,
            false);
    }
}

//=============================================================================
// External nested
//=============================================================================

void testExternalNested(
    STLGeometry& geometry)
{
    require(
        geometry.flowType ==
            FlowType::External,
        "external_nested requires external flow type.");


    bool hasChild =
        false;


    for(const auto& component :
        geometry.containment)
    {
        if(component.level == 1)
        {
            hasChild =
                true;

            break;
        }
    }


    require(
        hasChild,
        "external_nested requires at least one child component.");


    checkExternalSemantics(
        geometry);


    const OrientationSnapshot before =
        captureOrientation(
            geometry);

    const TopologySnapshot topologyBefore =
        captureTopology(
            geometry);


    geometry.interpretFlow();


    checkExternalSemantics(
        geometry);


    checkTopologyUnchanged(
        geometry,
        topologyBefore);


    for(std::size_t componentID = 0;
        componentID <
            geometry.topology.components.size();
        ++componentID)
    {
        checkComponentOrientationChange(
            geometry,
            before,
            componentID,
            false);
    }
}

//=============================================================================
// Internal to external
//=============================================================================

void testInternalToExternal(
    STLGeometry& geometry)
{
    require(
        geometry.flowType ==
            FlowType::Internal,
        "internal_to_external must start as internal.");


    require(
        countRoots(
            geometry) == 1,
        "internal_to_external requires exactly one root.");


    checkInternalSemantics(
        geometry);


    const OrientationSnapshot internalState =
        captureOrientation(
            geometry);

    const TopologySnapshot topologyBefore =
        captureTopology(
            geometry);


    geometry.flowType =
        FlowType::External;


    geometry.interpretFlow();


    checkExternalSemantics(
        geometry);


    checkTopologyUnchanged(
        geometry,
        topologyBefore);


    for(std::size_t componentID = 0;
        componentID <
            geometry.topology.components.size();
        ++componentID)
    {
        const std::size_t level =
            geometry.containment[
                componentID].level;


        if(level == 0)
        {
            checkComponentOrientationChange(
                geometry,
                internalState,
                componentID,
                true);
        }
        else
        {
            checkComponentOrientationChange(
                geometry,
                internalState,
                componentID,
                false);
        }
    }
}

//=============================================================================
// Test dispatch
//=============================================================================

void runTest(
    STLGeometry& geometry,
    const std::string& testCase)
{
    if(testCase ==
       "internal_root")
    {
        testInternalRoot(
            geometry);

        return;
    }


    if(testCase ==
       "internal_nested")
    {
        testInternalNested(
            geometry);

        return;
    }



    if(testCase ==
       "external_roots")
    {
        testExternalRoots(
            geometry);

        return;
    }


    if(testCase ==
       "external_nested")
    {
        testExternalNested(
            geometry);

        return;
    }


    if(testCase ==
       "internal_to_external")
    {
        testInternalToExternal(
            geometry);

        return;
    }


    throw std::runtime_error(
        "Unknown STL flow test case: " +
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


        if(testCase ==
           "internal_multiple_roots")
        {
            testInternalMultipleRoots(
                std::move(topology));
        }
        else
        {
            FlowType flowType =
                FlowType::Internal;


            if(testCase ==
                   "external_roots" ||
               testCase ==
                   "external_nested")
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
        }


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