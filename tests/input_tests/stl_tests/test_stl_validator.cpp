#include "stl_reader.hpp"
#include "stl_validator.hpp"

#include <exception>
#include <filesystem>
#include <iostream>
#include <string>

namespace {

using ntic::lbm::stl::STLData;


//=============================================================================
// Test helpers
//=============================================================================

bool expectValid(
    const std::filesystem::path& path,
    const std::string& testName)
{
    try {
        STLData data =
            ntic::lbm::stl::read(path);

        ntic::lbm::stl::validate(data);

        std::cout
            << "[PASS] "
            << testName
            << '\n';

        return true;
    }
    catch (const std::exception& e) {

        std::cerr
            << "[FAIL] "
            << testName
            << '\n'
            << "       Unexpected exception: "
            << e.what()
            << '\n';

        return false;
    }
}


bool expectInvalid(
    const std::filesystem::path& path,
    const std::string& testName)
{
    try {
        STLData data =
            ntic::lbm::stl::read(path);

        ntic::lbm::stl::validate(data);

        std::cerr
            << "[FAIL] "
            << testName
            << '\n'
            << "       Expected validation failure, "
            << "but validation succeeded."
            << '\n';

        return false;
    }
    catch (const std::exception& e) {

        std::cout
            << "[PASS] "
            << testName
            << '\n'
            << "       Detected: "
            << e.what()
            << '\n';

        return true;
    }
}

} // namespace


//=============================================================================
// Main
//=============================================================================

int main()
{
    const std::filesystem::path inputDir =
        TEST_INPUT_DIR;

    int failed = 0;

    //-------------------------------------------------------------------------
    // Valid facet
    //-------------------------------------------------------------------------

    if (!expectValid(
            inputDir / "valid_triangle.stl",
            "valid triangle")) {
        ++failed;
    }


    //-------------------------------------------------------------------------
    // Coincident vertices
    //-------------------------------------------------------------------------

    if (!expectInvalid(
            inputDir / "coincident_vertices.stl",
            "coincident vertices")) {
        ++failed;
    }


    //-------------------------------------------------------------------------
    // Nearly coincident vertices
    //-------------------------------------------------------------------------

    if (!expectInvalid(
            inputDir / "nearly_coincident_vertices.stl",
            "nearly coincident vertices")) {
        ++failed;
    }


    //-------------------------------------------------------------------------
    // Collinear vertices
    //-------------------------------------------------------------------------

    if (!expectInvalid(
            inputDir / "collinear_vertices.stl",
            "collinear vertices")) {
        ++failed;
    }


    //-------------------------------------------------------------------------
    // Nearly collinear vertices
    //-------------------------------------------------------------------------

    if (!expectInvalid(
            inputDir / "nearly_collinear_vertices.stl",
            "nearly collinear vertices")) {
        ++failed;
    }


    //-------------------------------------------------------------------------
    // Small but geometrically valid triangle
    //-------------------------------------------------------------------------

    if (!expectValid(
            inputDir / "small_valid_triangle.stl",
            "small valid triangle")) {
        ++failed;
    }


    //-------------------------------------------------------------------------
    // Summary
    //-------------------------------------------------------------------------

    std::cout << '\n';

    if (failed == 0) {
        std::cout
            << "All STL degenerate-facet tests passed."
            << '\n';

        return 0;
    }

    std::cerr
        << failed
        << " STL degenerate-facet test(s) failed."
        << '\n';

    return 1;
}