#include "circle.hpp"
#include "flow_type.hpp"
#include "vtk_output.hpp"
#include "vtk_output_2d.hpp"

#include <filesystem>
#include <iostream>
#include <stdexcept>


using namespace ntic::lbm::geometry;
using namespace ntic::lbm::geometry::test;


int main()
{
    try
    {
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
            << "========================================\n"
            << "Circle Geometry Test\n"
            << "========================================\n\n"

            << "Circle:\n"
            << "  center = ("
            << circle.center[0] << ", "
            << circle.center[1] << ")\n"
            << "  radius = "
            << circle.radius
            << "\n\n"

            << "Bounding box:\n"
            << "  min = ("
            << boundingBox.min[0] << ", "
            << boundingBox.min[1] << ")\n"
            << "  max = ("
            << boundingBox.max[0] << ", "
            << boundingBox.max[1] << ")\n\n"

            << "Open box:\n"
            << "  min = ("
            << openBox.min[0] << ", "
            << openBox.min[1] << ")\n"
            << "  max = ("
            << openBox.max[0] << ", "
            << openBox.max[1] << ")\n\n"

            << "Expected internal fluid region:\n"
            << "  Strictly inside the circle.\n"
            << "  Circle boundary is not fluid.\n\n"

            << "Expected external fluid region:\n"
            << "  Strictly inside the open box and\n"
            << "  strictly outside the circle.\n"
            << "  Circle boundary is not fluid.\n"
            << "  Open-box boundary is not fluid.\n\n";


        double spacing =
            0.0;

        std::cout
            << "Enter grid spacing: ";

        std::cin
            >> spacing;


        if(!std::cin ||
           spacing <= 0.0)
        {
            throw std::invalid_argument(
                "Grid spacing must be positive.");
        }


        const std::filesystem::path outputDirectory =
            "test_circle_vtks";

        const std::filesystem::path internalDirectory =
            outputDirectory /
            "internal_vtks";

        const std::filesystem::path externalDirectory =
            outputDirectory /
            "external_vtks";


        recreateOutputDirectory(
            outputDirectory);

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
                return
                    circle.contains(
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
                return
                    circle.contains(
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
            << "\nVTK output completed.\n\n"

            << "Internal files:\n"
            << "  "
            << (internalDirectory /
                "geometry.vtk")
            << "\n"
            << "  "
            << (internalDirectory /
                "bounding_box.vtk")
            << "\n\n"

            << "External files:\n"
            << "  "
            << (externalDirectory /
                "geometry.vtk")
            << "\n"
            << "  "
            << (externalDirectory /
                "bounding_box.vtk")
            << "\n"
            << "  "
            << (externalDirectory /
                "open_box.vtk")
            << "\n";


        return 0;
    }
    catch(const std::exception& error)
    {
        std::cerr
            << "Error: "
            << error.what()
            << "\n";

        return 1;
    }
}