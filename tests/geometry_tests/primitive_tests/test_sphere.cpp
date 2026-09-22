#include "sphere.hpp"
#include "vtk_output.hpp"
#include "vtk_output_3d.hpp"

#include <exception>
#include <filesystem>
#include <iostream>
#include <stdexcept>

using namespace ntic::lbm::geometry;
using namespace ntic::lbm::geometry::test;

int main()
{
    try
    {
        std::cout
            << "========================================\n"
            << "Sphere Geometry Test\n"
            << "========================================\n\n";

        const Sphere sphere{
            {4.0, 3.0, 5.0},
            2.0
        };

        const BoundingBox boundingBox =
            sphere.boundingBox();

        const BoundingBox openBox{
            {-1.0, -2.0, 0.0},
            {9.0, 8.0, 10.0}
        };

        std::cout
            << "Sphere:\n"
            << "  center = ("
            << sphere.center[0] << ", "
            << sphere.center[1] << ", "
            << sphere.center[2] << ")\n"
            << "  radius = "
            << sphere.radius
            << "\n\n";

        std::cout
            << "Bounding box:\n"
            << "  min = ("
            << boundingBox.min[0] << ", "
            << boundingBox.min[1] << ", "
            << boundingBox.min[2] << ")\n"
            << "  max = ("
            << boundingBox.max[0] << ", "
            << boundingBox.max[1] << ", "
            << boundingBox.max[2] << ")\n\n";

        std::cout
            << "Open box:\n"
            << "  min = ("
            << openBox.min[0] << ", "
            << openBox.min[1] << ", "
            << openBox.min[2] << ")\n"
            << "  max = ("
            << openBox.max[0] << ", "
            << openBox.max[1] << ", "
            << openBox.max[2] << ")\n\n";

        std::cout
            << "Expected internal fluid region:\n"
            << "  Strictly inside the sphere.\n"
            << "  Sphere surface is not fluid.\n\n";

        std::cout
            << "Expected external fluid region:\n"
            << "  Strictly inside the open box and\n"
            << "  strictly outside the sphere.\n"
            << "  Sphere surface is not fluid.\n"
            << "  Open-box surface is not fluid.\n\n";

        double spacing = 0.0;

        std::cout
            << "Enter grid spacing: ";

        std::cin >> spacing;

        if(!std::cin)
        {
            throw std::runtime_error(
                "Failed to read grid spacing.");
        }

        const std::filesystem::path outputDirectory =
            std::filesystem::path(
                GEOMETRY_TEST_OUTPUT_DIR) /
            "test_sphere_vtks";

        recreateOutputDirectory(
            outputDirectory);

        const std::filesystem::path internalDirectory =
            outputDirectory /
            "internal_vtks";

        const std::filesystem::path externalDirectory =
            outputDirectory /
            "