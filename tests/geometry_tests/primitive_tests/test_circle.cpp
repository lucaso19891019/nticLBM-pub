#include "circle.hpp"
#include "vtk_output.hpp"
#include "vtk_output_2d.hpp"

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
            << "Circle Geometry Test\n"
            << "========================================\n\n";

        const Circle circle{
            {4.0, 3.0, 0.0},
            2.0
        };

        const BoundingBox boundingBox =
            circle.boundingBox();

        const BoundingBox openBox{
            {-1.0, -2.0, 0.0},
            {9.0, 8.0, 0.0}
        };

        std::cout
            << "Circle:\n"
            << "  center = ("
            << circle.center[0] << ", "
            << circle.center[1] << ")\n"
            << "  radius = "
            << circle.radius
            << "\n\n";

        std::cout
            << "Bounding box:\n"
            << "  min = ("
            << boundingBox.min[0] << ", "
            << boundingBox.min[1] << ")\n"
            << "  max = ("
            << boundingBox.max[0] << ", "
            << boundingBox.max[1] << ")\n\n";

        std::cout
            << "Open box:\n"
            << "  min = ("
            << openBox.min[0] << ", "
            << openBox.min[1] << ")\n"
            << "  max = ("
            << openBox.max[0] << ", "
            << openBox.max[1] << ")\n\n";

        std::cout
            << "Expected internal fluid region:\n"
            << "  Strictly inside the circle.\n"
            << "  Circle boundary is not fluid.\n\n";

        std::cout
            << "Expected external fluid region:\n"
            << "  Strictly inside the open box and\n"
            << "  strictly outside the circle.\n"
            << "  Circle boundary is not fluid.\n"
            << "  Open-box boundary is not fluid.\n\n";

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
            "test_circle_vtks";

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

        writeVTK2DGeometry(
            internalDirectory /
                "geometry.vtk",
            boundingBox,
            spacing,
            [&circle](
                const Point& point)
            {
                return circle.contains(
                    point,
                    FlowType::Internal);
            });

        writeVTK2DBoundingBox(
            internalDirectory /
                "bounding_box.vtk",
            boundingBox);

        writeVTK2DGeometry(
            externalDirectory /
                "geometry.vtk",
            openBox,
            spacing,
            [&circle, &openBox](
                const Point& point)
            {
                return circle.contains(
                    point,
                    FlowType::External,
                    &openBox);
            });

        writeVTK2DBoundingBox(
            externalDirectory /
                "bounding_box.vtk",
            boundingBox);

        writeVTK2DBoundingBox(
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