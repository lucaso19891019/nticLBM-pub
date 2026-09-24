#include "circle.hpp"
#include "vtk_output.hpp"
#include "vtk_output_2d.hpp"

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
    const Circle& circle,
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


    const std::size_t pointCount =
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

        const std::size_t j =
            pointIndex /
            nx;

        const std::size_t i =
            pointIndex %
            nx;


        const Point point{
            domainBounds.min[0] +
                static_cast<double>(i) *
                spacing,

            domainBounds.min[1] +
                static_cast<double>(j) *
                spacing,

            domainBounds.min[2]
        };


        scalar[pointIndex] =
            circle.contains(
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


        const std::vector<double> internalScalar =
            buildScalar(
                circle,
                boundingBox,
                spacing,
                FlowType::Internal,
                nullptr);


        writeVTK2DScalar(
            internalDirectory /
                "geometry.vtk",
            internalNx,
            internalNy,
            boundingBox.min,
            spacing,
            internalScalar);


        writeVTK2DBoundingBox(
            internalDirectory /
                "bounding_box.vtk",
            boundingBox);


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


        const std::vector<double> externalScalar =
            buildScalar(
                circle,
                openBox,
                spacing,
                FlowType::External,
                &openBox);


        writeVTK2DScalar(
            externalDirectory /
                "geometry.vtk",
            externalNx,
            externalNy,
            openBox.min,
            spacing,
            externalScalar);


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