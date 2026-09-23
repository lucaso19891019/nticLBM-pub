#include "stl_geometry.hpp"

#include <cstddef>
#include <limits>
#include <stdexcept>
#include <utility>


namespace ntic::lbm::geometry
{

namespace
{

//=============================================================================
// Containment availability
//=============================================================================

void validateContainmentAvailable(
    const STLGeometry& geometry)
{
    const std::size_t componentCount =
        geometry.topology.components.size();


    if(componentCount == 0)
    {
        throw std::runtime_error(
            "Cannot interpret STL flow without components.");
    }


    if(geometry.containment.size() !=
       componentCount)
    {
        throw std::runtime_error(
            "STL component containment must be analyzed "
            "before flow interpretation.");
    }
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


    for(const STLComponentContainment& component :
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
// Flow structure validation
//=============================================================================

void validateFlowStructure(
    const FlowType flowType,
    const std::size_t rootCount)
{
    if(flowType ==
       FlowType::Internal)
    {
        if(rootCount != 1)
        {
            throw std::runtime_error(
                "Internal STL flow requires exactly "
                "one root component.");
        }


        return;
    }


    if(flowType ==
       FlowType::External)
    {
        if(rootCount == 0)
        {
            throw std::runtime_error(
                "External STL flow requires at least "
                "one root component.");
        }


        return;
    }


    throw std::runtime_error(
        "Unsupported STL flow type.");
}


//=============================================================================
// Internal flow semantics
//=============================================================================

void assignInternalFlow(
    STLGeometry& geometry)
{
    const std::size_t componentCount =
        geometry.containment.size();


    for(std::size_t componentID = 0;
        componentID < componentCount;
        ++componentID)
    {
        const std::size_t level =
            geometry.containment[
                componentID].level;


        STLComponentFlow& flow =
            geometry.flow[
                componentID];


        flow.active =
            true;


        if(level == 0)
        {
            flow.fluidSide =
                FluidSide::Inside;
        }
        else if(level == 1)
        {
            flow.fluidSide =
                FluidSide::Outside;
        }
        else
        {
            throw std::runtime_error(
                "Unsupported internal STL "
                "containment level.");
        }
    }
}


//=============================================================================
// External flow semantics
//=============================================================================

void assignExternalFlow(
    STLGeometry& geometry)
{
    const std::size_t componentCount =
        geometry.containment.size();


    for(std::size_t componentID = 0;
        componentID < componentCount;
        ++componentID)
    {
        const std::size_t level =
            geometry.containment[
                componentID].level;


        STLComponentFlow& flow =
            geometry.flow[
                componentID];


        if(level == 0)
        {
            flow.active =
                true;

            flow.fluidSide =
                FluidSide::Outside;
        }
        else if(level == 1)
        {
            flow.active =
                false;
        }
        else
        {
            throw std::runtime_error(
                "Unsupported external STL "
                "containment level.");
        }
    }
}


//=============================================================================
// Flow semantics
//=============================================================================

void assignFlowSemantics(
    STLGeometry& geometry,
    const FlowType flowType)
{
    if(flowType ==
       FlowType::Internal)
    {
        assignInternalFlow(
            geometry);

        return;
    }


    if(flowType ==
       FlowType::External)
    {
        assignExternalFlow(
            geometry);

        return;
    }


    throw std::runtime_error(
        "Unsupported STL flow type.");
}


//=============================================================================
// Component orientation
//=============================================================================

bool componentIsOutward(
    const STLGeometry& geometry,
    const std::size_t componentID)
{
    const double signedVolume =
        geometry.topology.components[
            componentID].signedVolume;


    if(signedVolume == 0.0)
    {
        throw std::runtime_error(
            "Cannot determine STL component orientation "
            "from zero signed volume.");
    }


    return
        signedVolume > 0.0;
}


//=============================================================================
// Flip component orientation
//=============================================================================

void flipComponent(
    STLGeometry& geometry,
    const std::size_t componentID)
{
    auto& topology =
        geometry.topology;

    auto& component =
        topology.components[
            componentID];


    const std::size_t facetCount =
        component.facets.size();


#pragma omp parallel for schedule(static)
    for(std::ptrdiff_t index = 0;
        index <
            static_cast<std::ptrdiff_t>(
                facetCount);
        ++index)
    {
        const std::size_t facetID =
            component.facets[
                static_cast<std::size_t>(
                    index)];


        auto& vertexIDs =
            topology.geometry.facetVertexIDs[
                facetID];


        std::swap(
            vertexIDs[1],
            vertexIDs[2]);


        auto& normal =
            topology.facetGeometry[
                facetID].normal;


        normal[0] =
            -normal[0];

        normal[1] =
            -normal[1];

        normal[2] =
            -normal[2];
    }


    component.signedVolume =
        -component.signedVolume;
}


//=============================================================================
// Orientation normalization
//=============================================================================

void normalizeOrientation(
    STLGeometry& geometry)
{
    const std::size_t componentCount =
        geometry.topology.components.size();


    for(std::size_t componentID = 0;
        componentID < componentCount;
        ++componentID)
    {
        const STLComponentFlow& flow =
            geometry.flow[
                componentID];


        if(!flow.active)
        {
            continue;
        }


        const bool currentlyOutward =
            componentIsOutward(
                geometry,
                componentID);


        const bool targetOutward =
            flow.fluidSide ==
                FluidSide::Outside;


        if(currentlyOutward !=
           targetOutward)
        {
            flipComponent(
                geometry,
                componentID);
        }
    }
}

} // namespace


//=============================================================================
// STL flow interpretation
//=============================================================================

void STLGeometry::interpretFlow(
    const FlowType flowType)
{
    validateContainmentAvailable(
        *this);


    const std::size_t rootCount =
        countRoots(
            *this);


    validateFlowStructure(
        flowType,
        rootCount);


    flow.clear();

    flow.resize(
        topology.components.size());


    assignFlowSemantics(
        *this,
        flowType);


    normalizeOrientation(
        *this);
}

} // namespace ntic::lbm::geometry