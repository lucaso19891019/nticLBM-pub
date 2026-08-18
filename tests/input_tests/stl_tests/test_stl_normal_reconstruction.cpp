#include "stl_reader.hpp"
#include "stl_validator.hpp"

#include <cmath>
#include <iostream>

namespace {

constexpr double TOLERANCE = 1.0e-12;

bool nearlyEqual(double a, double b)
{
    return std::abs(a - b) <= TOLERANCE;
}

} // namespace


int main(int argc, char* argv[])
{
    //-------------------------------------------------------------------------
    // Command-line arguments
    //-------------------------------------------------------------------------

    if (argc != 2) {
        std::cerr
            << "Usage: "
            << argv[0]
            << " <stl_file>\n";

        return 1;
    }


    //-------------------------------------------------------------------------
    // Read STL
    //-------------------------------------------------------------------------

    auto data =
        ntic::lbm::stl::read(argv[1]);

    if (data.facets.empty()) {
        std::cerr
            << "[FAIL] STL contains no facets.\n";

        return 1;
    }


    //-------------------------------------------------------------------------
    // Validate STL.
    //
    // validateFacetNormals() should replace the stored STL normal with the
    // unit geometric normal reconstructed from the vertex winding.
    //-------------------------------------------------------------------------

    try {
        ntic::lbm::stl::STLComponents components;
        ntic::lbm::stl::validate(data,"test",components);
    }
    catch (const std::exception& e) {
        std::cerr
            << "[FAIL] STL validation failed:\n"
            << "       "
            << e.what()
            << '\n';

        return 1;
    }


    //-------------------------------------------------------------------------
    // Check reconstructed normal.
    //
    // The test STL uses:
    //
    //     v0 = (0, 0, 0)
    //     v1 = (1, 0, 0)
    //     v2 = (0, 1, 0)
    //
    // Therefore:
    //
    //     (v1 - v0) x (v2 - v0) = (0, 0, 1)
    //
    // The original stored normal may have a different magnitude, but after
    // validation the final facet normal must be the unit vector (0, 0, 1).
    //-------------------------------------------------------------------------

    const auto& normal =
        data.facets[0].normal;

    if (!nearlyEqual(normal[0], 0.0) ||
        !nearlyEqual(normal[1], 0.0) ||
        !nearlyEqual(normal[2], 1.0)) {

        std::cerr
            << "[FAIL] Normal reconstruction failed.\n"
            << "       Expected: (0, 0, 1)\n"
            << "       Actual:   ("
            << normal[0] << ", "
            << normal[1] << ", "
            << normal[2] << ")\n";

        return 1;
    }


    //-------------------------------------------------------------------------
    // Test passed
    //-------------------------------------------------------------------------

    std::cout
        << "[PASS] STL facet normal reconstruction\n"
        << "       Normal: ("
        << normal[0] << ", "
        << normal[1] << ", "
        << normal[2] << ")\n";

    return 0;
}
