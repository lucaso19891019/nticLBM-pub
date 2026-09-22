#include "cylinder.hpp"
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
            << "Cylinder Geometry Test\n"
            << "========================================\n\n";

        const Cylinder cylinder{
            {4.0, 3.0, 5.0},
            2.0,
            6.0,
            Axis::Z
        };

        const BoundingBox boundingBox =
            cylinder.boundingBox();

        const BoundingBox openBox{
            {-1.0, -2.0, -1.0},
            {9.0, 8.0, 11.0}
        };

        std::cout
            << "Cylinder:\n"
            << "  center = ("
            << cylinder.center[0] << ", "
            << cylinder.center[1] << ", "
            << cylinder.center[2] << ")\n"
            << "  radius = "
            << cylinder.radius
            << "\n"
            << "  length = "
            << cylinder.length
            << "\n"
            << "  axis   = Z\n\n";

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
            << "  Strictly inside the cylinder.\n"
            << "  Cylinder surface is not fluid.\n\n";

        std::cout
            << "Expected external fluid region:\n"
            << "  Strictly inside the open box and\n"
            << "  strictly outside the cylinder.\n"
            << "  Cylinder surface is not fluid.\n"
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
            "test_cylinder_vtks";

        recreateOutputDirectory(
            outputDirectory);

        const std::filesystem::path internalDirectory =
            outputDirectory /
            "internal_vtks";

        const std::filesystem::path externalDirectory =
            outputDirectory /
            "external_vtks";

        std::filesystem::create_directories(
            internalDirectory);

        std::filesystem::create_directories(
            externalDirectory);

        writeVTK3DGeometry(
            internalDirectory /
                "geometry.vtk",
            boundingBox,
            spacing,
            [&cylinder](
                const Point& point)
            {
                return cylinder.contains(
                    point,
                    FlowType::Internal);
            });

        writeVTK3DBoundingBox(
            internalDirectory /
                "bounding_box.vtk",
            boundingBox);

        writeVTK3DGeometry(
            externalDirectory /
                "geometry.vtk",
            openBox,
            spacing,
            [&cylinder, &openBox](
                const Point& point)
            {
                return cylinder.contains(
                    point,
                    FlowType::External,
                    &openBox);
            });

        writeVTK3DBoundingBox(
            externalDirectory /
                "bounding_box.vtk",
            boundingBox);

        writeVTK3DBoundingBox(
            externalDirectory /
                "open_box.vtk",
            openBox);

        std::cout
            << "\nVTK output written to:\n"
            << outputDirectory
            << "\n";

        return 0;
    }
    catch(const std::exception& exception)
    {
        std::cerr
            << "Error: "
            << exception.what()
            << "\n";

        return 1;
    }
}