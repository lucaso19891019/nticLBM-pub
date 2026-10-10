#include "stl_geometry.hpp"

#include <cstddef>
#include <iostream>
#include <stdexcept>
#include <utility>


namespace ntic::lbm::geometry
{

namespace
{

//=============================================================================
// STL bounding box
//=============================================================================

BoundingBox makeSTLBounds(
    const STLGeometry& geometry)
{
    const std::size_t vertexCount =
        geometry.topology.geometry.vertices.size();


    BoundingBox bounds;


#pragma omp parallel
    {
        BoundingBox localBounds;

        bool hasVertex =
            false;


#pragma omp for schedule(static)
        for(std::ptrdiff_t index = 0;
            index <
                static_cast<std::ptrdiff_t>(
                    vertexCount);
            ++index)
        {
            localBounds.expand(
                geometry.topology.geometry.vertices[
                    static_cast<std::size_t>(
                        index)]);

            hasVertex =
                true;
        }


        if(hasVertex)
        {
#pragma omp critical
            {
                bounds.expand(
                    localBounds.min);

                bounds.expand(
                    localBounds.max);
            }
        }
    }


    return bounds;
}


//=============================================================================
// Open-box validation
//=============================================================================

void validateOpenBox(
    const BoundingBox& openBox,
    const BoundingBox& bounds)
{
    if(openBox.max[0] <= openBox.min[0] ||
       openBox.max[1] <= openBox.min[1] ||
       openBox.max[2] <= openBox.min[2])
    {
        throw std::runtime_error(
            "External STL open box must have positive dimensions.");
    }


    if(openBox.min[0] > bounds.min[0] ||
       openBox.min[1] > bounds.min[1] ||
       openBox.min[2] > bounds.min[2] ||
       openBox.max[0] < bounds.max[0] ||
       openBox.max[1] < bounds.max[1] ||
       openBox.max[2] < bounds.max[2])
    {
        throw std::runtime_error(
            "External STL open box must contain the complete STL bounding box.");
    }
}


//=============================================================================
// Interactive open-box input
//=============================================================================

BoundingBox readOpenBox(
    const BoundingBox& bounds)
{
    std::cout
        << "External flow requires an open box.\n\n"
        << "STL tight bounding box:\n"
        << "  min    = ("
        << bounds.min[0] << ", "
        << bounds.min[1] << ", "
        << bounds.min[2] << ")\n"
        << "  max    = ("
        << bounds.max[0] << ", "
        << bounds.max[1] << ", "
        << bounds.max[2] << ")\n"
        << "  width  = "
        << bounds.width()
        << "\n"
        << "  height = "
        << bounds.height()
        << "\n"
        << "  depth  = "
        << bounds.depth()
        << "\n\n";


    Point origin;

    Point size;


    std::cout
        << "Enter open-box origin (x y z): ";

    std::cin
        >> origin[0]
        >> origin[1]
        >> origin[2];


    if(!std::cin)
    {
        throw std::runtime_error(
            "Failed to read external open-box origin.");
    }


    std::cout
        << "Enter open-box size (width height depth): ";

    std::cin
        >> size[0]
        >> size[1]
        >> size[2];


    if(!std::cin)
    {
        throw std::runtime_error(
            "Failed to read external open-box size.");
    }


    BoundingBox openBox;

    openBox.min =
        origin;

    openBox.max =
    {
        origin[0] + size[0],
        origin[1] + size[1],
        origin[2] + size[2]
    };


    validateOpenBox(
        openBox,
        bounds);


    return openBox;
}

} // namespace


//=============================================================================
// STL geometry construction
//=============================================================================

STLGeometry::STLGeometry(
    stl::FacetTopology inputTopology,
    const FlowType inputFlowType,
    const BoundingBox* inputOpenBox)
    : topology(
          std::move(
              inputTopology)),
      flowType(
          inputFlowType)
{
    if(topology.geometry.vertices.empty())
    {
        throw std::runtime_error(
            "Cannot construct empty STL geometry.");
    }


    bounds =
        makeSTLBounds(
            *this);


    if(flowType ==
       FlowType::External)
    {
        if(inputOpenBox != nullptr)
        {
            validateOpenBox(
                *inputOpenBox,
                bounds);

            openBox =
                *inputOpenBox;
        }
        else
        {
            openBox =
                readOpenBox(
                    bounds);
        }
    }


    analyzeContainment();
    
    interpretFlow();
    
    if(flowType == FlowType::Internal)
    {
        identifyBoundaryFeatures();
    }
}

} // namespace ntic::lbm::geometry
