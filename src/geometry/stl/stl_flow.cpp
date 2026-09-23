#include "stl_geometry.hpp"

#include <cstddef>
#include <limits>
#include <stdexcept>
#include <utility>


namespace ntic::lbm::geometry
{

//=============================================================================
// Containment availability
//=============================================================================

void STLGeometry::validateContainmentAvailable() const
{
    const std::size_t componentCount =
        topology.components.size();


    if(componentCount == 0)
    {
        throw std::runtime_error(
            "Cannot interpret STL flow without components.");
    }


    if(containment.size() !=
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

std::size_t STLGeometry::countRoots() const
{
    const std::size_t noParent =
        std::numeric_limits<std::size_t>::max();


    std::size_t rootCount =
        0;


    for(const STLComponentContainment& component :
        containment)
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

void STLGeometry::validateFlowStructure(
    const std::size_t rootCount) const
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

void STLGeometry::assignInternalFlow()
{
    const std::size_t componentCount =
        containment.size();


    for(std::size_t componentID = 0;
        componentID < componentCount;
        ++componentID)
    {
        const std::size_t level =
            containment[
                componentID].level;


        STLComponentFlow& componentFlow =
            flow[
                componentID];


        componentFlow.active =
            true;


        if(level == 0)
        {
            componentFlow.fluidSide =
                FluidSide::Inside;
        }
        else if(level == 1)
        {
            componentFlow.fluidSide =
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

void STLGeometry::assignExternalFlow()
{
    const std::size_t componentCount =
        containment.size();


    for(std::size_t componentID = 0;
        componentID < componentCount;
        ++componentID)
    {
        const std::size_t level =
            containment[
                componentID].level;


        STLComponentFlow& componentFlow =
            flow[
                componentID];


        if(level == 0)
        {
            componentFlow.active =
                true;

            componentFlow.fluidSide =
                FluidSide::Outside;
        }
        else if(level == 1)
        {
            componentFlow.active =
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

void STLGeometry::assignFlowSemantics()
{
    if(flowType ==
       FlowType::Internal)
    {
        assignInternalFlow();

        return;
    }


    if(flowType ==
       FlowType::External)
    {
        assignExternalFlow();

        return;
    }


    throw std::runtime_error(
        "Unsupported STL flow type.");
}


//=============================================================================
// Component orientation
//=============================================================================

bool STLGeometry::componentIsOutward(
    const std::size_t componentID) const
{
    const double signedVolume =
        topology.components[
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

void STLGeometry::flipComponent(
    const std::size_t componentID)
{
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

void STLGeometry::normalizeOrientation()
{
    const std::size_t componentCount =
        topology.components.size();


    for(std::size_t componentID = 0;
        componentID < componentCount;
        ++componentID)
    {
        const STLComponentFlow& componentFlow =
            flow[
                componentID];


        if(!componentFlow.active)
        {
            continue;
        }


        const bool currentlyOutward =
            componentIsOutward(
                componentID);


        const bool targetOutward =
            componentFlow.fluidSide ==
                FluidSide::Outside;


        if(currentlyOutward !=
           targetOutward)
        {
            flipComponent(
                componentID);
        }
    }
}


//=============================================================================
// STL flow interpretation
//=============================================================================

void STLGeometry::interpretFlow()
{
    validateContainmentAvailable();


    const std::size_t rootCount =
        countRoots();


    validateFlowStructure(
        rootCount);


    flow.clear();

    flow.resize(
        topology.components.size());


    assignFlowSemantics();


    normalizeOrientation();
}

} // namespace ntic::lbm::geometry