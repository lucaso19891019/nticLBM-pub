#include "box.hpp"
#include "flow_type.hpp"
#include "vtk_output.hpp"
#include "vtk_output_3d.hpp"

#include <filesystem>
#include <iostream>
#include <stdexcept>


using namespace ntic::lbm::geometry;
using namespace ntic::lbm::geometry::test;


int main()
{
    try
    {
        const Box box{
            {1.0, 2.0, 3.0},
            {7.0, 6.0, 8.0}
        };


        const BoundingBox boundingBox =
            box.boundingBox();


        const BoundingBox openBox{
            {-2.0, -1.0, 0.0},
            {10.0, 9.0, 11.0}
        };


        std::cout
            << "========================================\n"
            << "Box Geometry Test\n"
            << "========================================\n\n"

            << "Box:\n"
            << "  min    = ("
            << box.min[0] << ", "
            << box.min[1] << ", "
            << box.min[2] << ")\n"
            << "  max    = ("
            << box.max[0] << ", "
            << box.max[1] << ", "
            << box.max[2] << ")\n"
            << "  width  = "
            << box.max[0] -
               box.min[0]
            << "\n"
            << "  height = "
            << box.max[1] -
               box.min[1]
            << "\n"
            << "  depth  = "
            << box.max[2] -
               box.min[2]
            << "\n\n"

            << "Bounding box:\n"
            << "  min = ("
            << boundingBox.min[0] << ", "
            << boundingBox.min[1] << ", "
            << boundingBox.min[2] << ")\n"
            << "  max = ("
            << boundingBox.max[0] << ", "
            << boundingBox.max[1] << ", "
            << boundingBox.max[2] << ")\n\n"

            << "Open box:\n"
            << "  min = ("
            << openBox.min[0] << ", "
            << openBox.min[1] << ", "
            << openBox.min[2] << ")\n"
            << "  max = ("
            << openBox.max[0] << ", "
            << openBox.max[1] << ", "
            << openBox.max[2] << ")\n\n"

            << "Expected internal fluid region:\n"
            << "  Strictly inside the box.\n"
            << "  Box surface is not fluid.\n\n"

            << "Expected external fluid region:\n"
            << "  Strictly inside the open box and\n"
            << "  strictly outside the box.\n"
            << "  Box surface is not fluid.\n"
            << "  Open-box surface is not fluid.\n\n";


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
            "test_box_vtks";

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


        writeVTK3DGeometry(
            internalDirectory /
                "geometry.vtk",
            boundingBox,
            spacing,
            [&box](
                const Point& point)
            {
                return
                    box.contains(
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
            [&box, &openBox](
                const Point& point)
            {
                return
                    box.contains(
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