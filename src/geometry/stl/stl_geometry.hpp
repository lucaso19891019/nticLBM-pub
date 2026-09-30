#pragma once

#include "bounding_box.hpp"
#include "cell_type.hpp"
#include "geometry_analysis.hpp"
#include "flow_type.hpp"
#include "point.hpp"
#include "stl_topology.hpp"

#include <cstddef>
#include <limits>
#include <vector>


namespace ntic::lbm::geometry
{

//=============================================================================
// STL component containment
//=============================================================================

struct STLComponentContainment
{
    std::size_t parent =
        std::numeric_limits<std::size_t>::max();

    std::vector<std::size_t> children;

    std::size_t root =
        std::numeric_limits<std::size_t>::max();

    std::size_t level =
        0;
};

//=============================================================================
// STL component flow
//=============================================================================

enum class FluidSide
{
    Inside,
    Outside
};


struct STLComponentFlow
{
    bool active =
        false;

    FluidSide fluidSide =
        FluidSide::Outside;
};


//=============================================================================
// STLGeometry
//=============================================================================

struct STLGeometry
{
    stl::FacetTopology topology;

    BoundingBox bounds;

    BoundingBox openBox;

    std::vector<STLComponentContainment> containment;

    std::vector<STLComponentFlow> flow;

    FlowType flowType;


    STLGeometry(
        stl::FacetTopology topology,
        FlowType flowType = FlowType::Internal,
        const BoundingBox* openBox = nullptr);


    //-------------------------------------------------------------------------
    // Containment
    //-------------------------------------------------------------------------

    void analyzeContainment();


    [[nodiscard]]
    bool pointInComponent(
        std::size_t componentID,
        const Point& point) const;


    [[nodiscard]]
    Point componentTestPoint(
        std::size_t componentID) const;


    [[nodiscard]]
    bool componentInComponent(
        std::size_t innerComponentID,
        std::size_t outerComponentID) const;


    //-------------------------------------------------------------------------
    // Flow
    //-------------------------------------------------------------------------

    void interpretFlow();


    void validateContainmentAvailable() const;


    [[nodiscard]]
    std::size_t countRoots() const;


    void validateFlowStructure(
        std::size_t rootCount) const;


    void assignInternalFlow();


    void assignExternalFlow();


    void assignFlowSemantics();


    [[nodiscard]]
    bool componentIsOutward(
        std::size_t componentID) const;


    void flipComponent(
        std::size_t componentID);


    void normalizeOrientation();


    //-------------------------------------------------------------------------
    // Translation
    //-------------------------------------------------------------------------

    void translate(
        const Point* targetPoint = nullptr);


    //-------------------------------------------------------------------------
    // Interior area analysis
    //-------------------------------------------------------------------------

    template <typename LatticeModel>
    void analysis(
        GeometryAnalysis3D<LatticeModel>& analysis) const
    {
        interiorAreaAnalysis(
            analysis.gridSpacing,
            analysis.domain,
            analysis.nx,
            analysis.ny,
            analysis.nz,
            analysis.scalar,
            analysis.boundaryX,
            analysis.boundaryY,
            analysis.boundaryZ);

        boundaryQAnalysis(
            analysis);
    }


    template <typename LatticeModel>
    void boundaryQAnalysis(
        GeometryAnalysis3D<LatticeModel>&) const
    {
        // Reserved for lattice-link q analysis.
    }


    void interiorAreaAnalysis(
        double gridSpacing,
        BoundingBox& domain,
        std::size_t& nx,
        std::size_t& ny,
        std::size_t& nz,
        std::vector<double>& scalar,
        std::vector<double>& boundaryX,
        std::vector<double>& boundaryY,
        std::vector<double>& boundaryZ) const;
};

} // namespace ntic::lbm::geometry