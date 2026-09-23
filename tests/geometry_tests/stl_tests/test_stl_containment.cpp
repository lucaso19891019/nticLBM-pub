#include "flow_type.hpp"
#include "stl_geometry.hpp"
#include "stl_reader.hpp"
#include "stl_validator.hpp"

#include <cstddef>
#include <exception>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>


namespace
{

using ntic::lbm::geometry::FlowType;
using ntic::lbm::geometry::STLGeometry;
using ntic::lbm::geometry::STLComponentContainment;


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


void requireComponentCount(
    const STLGeometry& geometry,
    const std::size_t expected)
{
    require(
        geometry.containment.size() ==
            expected,
        "Unexpected containment component count.");
}


void requireRoot(
    const STLGeometry& geometry,
    const std::size_t componentID)
{
    const std::size_t noParent =
        std::numeric_limits<std::size_t>::max();


    const STLComponentContainment& component =
        geometry.containment[
            componentID];


    require(
        component.parent ==
            noParent,
        "Expected component to be a root.");

    require(
        component.root ==
            componentID,
        "Root component has incorrect root ID.");

    require(
        component.level ==
            0,
        "Root component has incorrect level.");
}


void requireChild(
    const STLGeometry& geometry,
    const std::size_t componentID,
    const std::size_t parentID,
    const std::size_t rootID)
{
    const STLComponentContainment& component =
        geometry.containment[
            componentID];


    require(
        component.parent ==
            parentID,
        "Child component has incorrect parent.");

    require(
        component.root ==
            rootID,
        "Child component has incorrect root.");

    require(
        component.level ==
            1,
        "Child component has incorrect level.");
}


bool hasChild(
    const STLComponentContainment& component,
    const std::size_t childID)
{
    for(const std::size_t current :
        component.children)
    {
        if(current ==
           childID)
        {
            return true;
        }
    }


    return false;
}


//=============================================================================
// Single root
//=============================================================================

void testSingleRoot(
    const STLGeometry& geometry)
{
    requireComponentCount(
        geometry,
        1);


    requireRoot(
        geometry,
        0);


    require(
        geometry.containment[0].children.empty(),
        "Single root unexpectedly contains children.");
}


//=============================================================================
// Two roots
//=============================================================================

void testTwoRoots(
    const STLGeometry& geometry)
{
    requireComponentCount(
        geometry,
        2);


    requireRoot(
        geometry,
        0);

    requireRoot(
        geometry,
        1);


    require(
        geometry.containment[0].children.empty(),
        "First root unexpectedly contains children.");

    require(
        geometry.containment[1].children.empty(),
        "Second root unexpectedly contains children.");
}


//=============================================================================
// One root with one child
//=============================================================================

void testNested(
    const STLGeometry& geometry)
{
    requireComponentCount(
        geometry,
        2);


    std::size_t rootID =
        std::numeric_limits<std::size_t>::max();

    std::size_t childID =
        std::numeric_limits<std::size_t>::max();


    for(std::size_t componentID = 0;
        componentID < 2;
        ++componentID)
    {
        if(geometry.containment[
               componentID].level == 0)
        {
            rootID =
                componentID;
        }
        else if(geometry.containment[
                    componentID].level == 1)
        {
            childID =
                componentID;
        }
    }


    require(
        rootID !=
            std::numeric_limits<std::size_t>::max(),
        "Nested geometry contains no root.");

    require(
        childID !=
            std::numeric_limits<std::size_t>::max(),
        "Nested geometry contains no child.");


    requireRoot(
        geometry,
        rootID);

    requireChild(
        geometry,
        childID,
        rootID,
        rootID);


    require(
        geometry.containment[
            rootID].children.size() == 1,
        "Nested root must contain exactly one child.");

    require(
        hasChild(
            geometry.containment[
                rootID],
            childID),
        "Nested root does not reference its child.");

    require(
        geometry.containment[
            childID].children.empty(),
        "Nested child unexpectedly contains children.");
}


//=============================================================================
// One root with two children
//=============================================================================

void testTwoChildren(
    const STLGeometry& geometry)
{
    requireComponentCount(
        geometry,
        3);


    std::size_t rootID =
        std::numeric_limits<std::size_t>::max();


    for(std::size_t componentID = 0;
        componentID < 3;
        ++componentID)
    {
        if(geometry.containment[
               componentID].level == 0)
        {
            require(
                rootID ==
                    std::numeric_limits<std::size_t>::max(),
                "Expected exactly one root.");

            rootID =
                componentID;
        }
    }


    require(
        rootID !=
            std::numeric_limits<std::size_t>::max(),
        "Geometry contains no root.");


    requireRoot(
        geometry,
        rootID);


    require(
        geometry.containment[
            rootID].children.size() == 2,
        "Root must contain exactly two children.");


    std::size_t childCount =
        0;


    for(std::size_t componentID = 0;
        componentID < 3;
        ++componentID)
    {
        if(componentID ==
           rootID)
        {
            continue;
        }


        requireChild(
            geometry,
            componentID,
            rootID,
            rootID);


        require(
            hasChild(
                geometry.containment[
                    rootID],
                componentID),
            "Root does not reference expected child.");


        require(
            geometry.containment[
                componentID].children.empty(),
            "Sibling component unexpectedly contains children.");


        ++childCount;
    }


    require(
        childCount == 2,
        "Expected exactly two child components.");
}


//=============================================================================
// Two independent nested trees
//=============================================================================

void testTwoNestedRoots(
    const STLGeometry& geometry)
{
    requireComponentCount(
        geometry,
        4);


    std::size_t rootCount =
        0;

    std::size_t childCount =
        0;


    for(std::size_t componentID = 0;
        componentID < 4;
        ++componentID)
    {
        const STLComponentContainment& component =
            geometry.containment[
                componentID];


        if(component.level == 0)
        {
            ++rootCount;


            requireRoot(
                geometry,
                componentID);


            require(
                component.children.size() == 1,
                "Each root must contain exactly one child.");
        }
        else if(component.level == 1)
        {
            ++childCount;


            require(
                component.parent !=
                    std::numeric_limits<std::size_t>::max(),
                "Child has no parent.");

            require(
                component.root ==
                    component.parent,
                "Child root does not match its parent.");

            require(
                hasChild(
                    geometry.containment[
                        component.parent],
                    componentID),
                "Parent does not reference its child.");

            require(
                component.children.empty(),
                "Child unexpectedly contains children.");
        }
        else
        {
            throw std::runtime_error(
                "Unexpected containment level.");
        }
    }


    require(
        rootCount == 2,
        "Expected exactly two roots.");

    require(
        childCount == 2,
        "Expected exactly two children.");
}


//=============================================================================
// Expected three-level failure
//=============================================================================

void testThreeLevels(
    ntic::lbm::stl::FacetTopology topology)
{
    bool caughtExpectedException =
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
        caughtExpectedException =
            true;
    }


    require(
        caughtExpectedException,
        "Three-level containment did not throw.");
}


//=============================================================================
// Run successful containment case
//=============================================================================

void runSuccessfulCase(
    STLGeometry& geometry,
    const std::string& caseName)
{
    if(caseName ==
       "single_root")
    {
        testSingleRoot(
            geometry);
    }
    else if(caseName ==
            "two_roots")
    {
        testTwoRoots(
            geometry);
    }
    else if(caseName ==
            "aabb_overlap")
    {
        testTwoRoots(
            geometry);
    }
    else if(caseName ==
            "pseudo_containment")
    {
        testTwoRoots(
            geometry);
    }
    else if(caseName ==
            "nested")
    {
        testNested(
            geometry);
    }
    else if(caseName ==
            "two_children")
    {
        testTwoChildren(
            geometry);
    }
    else if(caseName ==
            "two_nested_roots")
    {
        testTwoNestedRoots(
            geometry);
    }
    else
    {
        throw std::runtime_error(
            "Unknown STL containment test case.");
    }
}


//=============================================================================
// Main
//=============================================================================

} // namespace


int main(
    const int argc,
    char* argv[])
{
    if(argc != 3)
    {
        std::cerr
            << "Usage: "
            << argv[0]
            << " <case> <stl_file>"
            << std::endl;

        return 1;
    }


    try
    {
        const std::string caseName =
            argv[1];

        const std::string stlFile =
            argv[2];


        ntic::lbm::stl::STLData data =
            ntic::lbm::stl::read(
                stlFile);


        ntic::lbm::stl::FacetTopology topology;


        ntic::lbm::stl::validate(
            data,
            "full",
            topology);


        if(caseName ==
           "three_levels")
        {
            testThreeLevels(
                std::move(topology));
        }
        else
        {
            FlowType flowType =
                FlowType::Internal;


            if(caseName ==
                   "two_roots" ||
               caseName ==
                   "aabb_overlap" ||
               caseName ==
                   "pseudo_containment" ||
               caseName ==
                   "two_nested_roots")
            {
                flowType =
                    FlowType::External;
            }


            STLGeometry geometry(
                std::move(topology),
                flowType);


            runSuccessfulCase(
                geometry,
                caseName);
        }


        std::cout
            << "[PASS] "
            << caseName
            << " : "
            << stlFile
            << std::endl;


        return 0;
    }
    catch(const std::exception& exception)
    {
        std::cerr
            << "[FAIL] "
            << exception.what()
            << std::endl;


        return 1;
    }
}