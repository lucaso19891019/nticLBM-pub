#include "rectangle.hpp"
#include "vtk_output.hpp"
#include "vtk_output_2d.hpp"

#include <exception>
#include <filesystem>
#include <iostream>

using namespace ntic::lbm::geometry;
using namespace ntic::lbm::geometry::test;

int main()
{
    try
    {
        std::cout
            << "========================================\n"
            << "Rectangle Geometry Test\n"
            << "========================================\n\n";

        const Rectangle rectangle{
            {1.0, 2.0, 0.0},
            {7.0, 6.0, 0.0}
        };

        const BoundingBox boundingBox =
            rectangle.boundingBox();

        const BoundingBox openBox{
            {-2.0, -1.0, 0.0},
            {10.0, 9.0, 0.0}
        };

        std::cout
            << "Rectangle:\n"
            << "  min = ("
            << rectangle.min[0] << ", "
            << rectangle.min[1] << ")\n"
            << "  max = ("
            << rectangle.max[0] << ", "
            << rectangle.max[1] << ")\n\n";

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
            << "  Strictly inside the rectangle.\n"
            << "  Rectangle boundary is not fluid.\n\n";

        std::cout
            << "Expected external fluid region:\n"
            << "  Strictly inside the open box and\n"
            << "  strictly outside the rectangle.\n"
            << "  Rectangle boundary is not fluid.\n"
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
            "test_rectangle_vtks";

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
            [&rectangle](
                const Point& point)
            {
                return rectangle.contains(
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
            [&rectangle, &openBox](
                const Point& point)
            {
                return rectangle.contains(
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