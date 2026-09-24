#include "cylinder.hpp"
#include "vtk_output.hpp"
#include "vtk_output_3d.hpp"

#include <cmath>
#include <cstddef>
#include <exception>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <vector>


using namespace ntic::lbm::geometry;
using namespace ntic::lbm::geometry::test;


namespace
{

std::vector<double> buildScalar(
    const Cylinder& cylinder,
    const BoundingBox& domainBounds,
    const double spacing,
    const FlowType flowType,
    const BoundingBox* openBox)
{
    const std::size_t nx =
        static_cast<std::size_t>(
            std::floor(
                domainBounds.width() /
                spacing)) +
        1;

    const std::size_t ny =
        static_cast<std::size_t>(
            std::floor(
                domainBounds.height() /
                spacing)) +
        1;

    const std::size_t nz =
        static_cast<std::size_t>(
            std::floor(
                domainBounds.depth() /
                spacing)) +
        1;


    const std::size_t pointCount =
        nx * ny * nz;

    const std::size_t xySize =
        nx * ny;


    std::vector<double> scalar(
        pointCount,
        0.0);


#pragma omp parallel for schedule(static)
    for(std::ptrdiff_t index = 0;
        index <
            static_cast<std::ptrdiff_t>(
                pointCount);
        ++index)
    {
        const std::size_t pointIndex =
            static_cast<std::size_t>(
                index);


        const std::size_t k =
            pointIndex /
            xySize;

        const std::size_t remainder =
            pointIndex %
            xySize;

        const std::size_t j =
            remainder /
            nx;

        const std::size_t i =
            remainder %
            nx;


        const Point point{
            domainBounds.min[0] +
                static_cast<double>(i) *
                spacing,

            domainBounds.min[1] +
                static_cast<double>(j) *
                spacing,

            domainBounds.min[2] +
                static_cast<double>(k) *
                spacing
        };


        scalar[pointIndex] =
            cylinder.contains(
                point,
                flowType,
                openBox)
                ? 1.0
                : 0.0;
    }


    return scalar;
}

} // namespace


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


        double spacing =
            0.0;


        std::cout
            << "Enter grid spacing: ";

        std::cin
            >> spacing;


        if(!std::cin)
        {
            throw std::runtime_error(
                "Failed to read grid spacing.");
        }


        if(spacing <= 0.0)
        {
            throw std::runtime_error(
                "Grid spacing must be positive.");
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


        //---------------------------------------------------------------------
        // Internal
        //---------------------------------------------------------------------

        const std::size_t internalNx =
            static_cast<std::size_t>(
                std::floor(
                    boundingBox.width() /
                    spacing)) +
            1;

        const std::size_t internalNy =
            static_cast<std::size_t>(
                std::floor(
                    boundingBox.height() /
                    spacing)) +
            1;

        const std::size_t internalNz =
            static_cast<std::size_t>(
                std::floor(
                    boundingBox.depth() /
                    spacing)) +
            1;


        const std::vector<double> internalScalar =
            buildScalar(
                cylinder,
                boundingBox,
                spacing,
                FlowType::Internal,
                nullptr);


        writeVTK3DScalar(
            internalDirectory /
                "geometry.vtk",
            internalNx,
            internalNy,
            internalNz,
            boundingBox.min,
            spacing,
            internalScalar);


        writeVTK3DBoundingBox(
            internalDirectory /
                "bounding_box.vtk",
            boundingBox);


        //---------------------------------------------------------------------
        // External
        //---------------------------------------------------------------------

        const std::size_t externalNx =
            static_cast<std::size_t>(
                std::floor(
                    openBox.width() /
                    spacing)) +
            1;

        const std::size_t externalNy =
            static_cast<std::size_t>(
                std::floor(
                    openBox.height() /
                    spacing)) +
            1;

        const std::size_t externalNz =
            static_cast<std::size_t>(
                std::floor(
                    openBox.depth() /
                    spacing)) +
            1;


        const std::vector<double> externalScalar =
            buildScalar(
                cylinder,
                openBox,
                spacing,
                FlowType::External,
                &openBox);


        writeVTK3DScalar(
            externalDirectory /
                "geometry.vtk",
            externalNx,
            externalNy,
            externalNz,
            openBox.min,
            spacing,
            externalScalar);


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