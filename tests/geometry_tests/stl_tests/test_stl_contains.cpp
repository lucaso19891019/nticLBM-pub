#include "stl_geometry.hpp"
#include "stl_reader.hpp"
#include "stl_validator.hpp"
#include "vtk_output.hpp"
#include "vtk_output_3d.hpp"

#include <cmath>
#include <cstddef>
#include <exception>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <utility>
#include <vector>


using namespace ntic::lbm::geometry;
using namespace ntic::lbm::geometry::test;


namespace
{

double cellTypeScalar(
    const CellType cellType)
{
    if(cellType ==
       CellType::Dry)
    {
        return 0.0;
    }


    if(cellType ==
       CellType::Boundary)
    {
        return 0.5;
    }


    return 1.0;
}


BoundingBox makeVisualizationBox(
    const BoundingBox& bounds,
    const double spacing)
{
    const double padding =
        2.0 *
        spacing;


    return
    {
        {
            bounds.min[0] -
                padding,
            bounds.min[1] -
                padding,
            bounds.min[2] -
                padding
        },
        {
            bounds.max[0] +
                padding,
            bounds.max[1] +
                padding,
            bounds.max[2] +
                padding
        }
    };
}


std::vector<double> buildScalar(
    const STLGeometry& geometry,
    const BoundingBox& domainBounds,
    const double spacing)
{
    const std::size_t nx =
        static_cast<std::size_t>(
            std::ceil(
                domainBounds.width() /
                spacing));

    const std::size_t ny =
        static_cast<std::size_t>(
            std::ceil(
                domainBounds.height() /
                spacing));

    const std::size_t nz =
        static_cast<std::size_t>(
            std::ceil(
                domainBounds.depth() /
                spacing));


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
                (
                    static_cast<double>(i) +
                    0.5
                ) *
                spacing,

            domainBounds.min[1] +
                (
                    static_cast<double>(j) +
                    0.5
                ) *
                spacing,

            domainBounds.min[2] +
                (
                    static_cast<double>(k) +
                    0.5
                ) *
                spacing
        };


        const CellType cellType =
            geometry.contains(
                point,
                spacing);


        scalar[pointIndex] =
            cellTypeScalar(
                cellType);
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
            << "STL Contains Test\n"
            << "========================================\n\n";


        const std::filesystem::path stlFile =
            "../tests/geometry_tests/stl_tests/stls/"
            "smooth_irregular_branched_channel.stl";

        //---------------------------------------------------------------------
        // Internal geometry
        //---------------------------------------------------------------------

        ntic::lbm::stl::STLData internalData =
            ntic::lbm::stl::read(
                stlFile.string());


        ntic::lbm::stl::FacetTopology internalTopology;


        ntic::lbm::stl::validate(
            internalData,
            "full",
            internalTopology);


        STLGeometry internalGeometry(
            std::move(
                internalTopology),
            FlowType::Internal);


        //---------------------------------------------------------------------
        // External geometry
        //---------------------------------------------------------------------

        ntic::lbm::stl::STLData externalData =
            ntic::lbm::stl::read(
                stlFile.string());


        ntic::lbm::stl::FacetTopology externalTopology;


        ntic::lbm::stl::validate(
            externalData,
            "full",
            externalTopology);


        STLGeometry externalGeometry(
            std::move(
                externalTopology),
            FlowType::External);

        //---------------------------------------------------------------------
        // STL geometry scale
        //---------------------------------------------------------------------

        std::cout
            << "STL bounding box:\n"
            << "  min    = ("
            << internalGeometry.bounds.min[0] << ", "
            << internalGeometry.bounds.min[1] << ", "
            << internalGeometry.bounds.min[2] << ")\n"
            << "  max    = ("
            << internalGeometry.bounds.max[0] << ", "
            << internalGeometry.bounds.max[1] << ", "
            << internalGeometry.bounds.max[2] << ")\n"
            << "  width  = "
            << internalGeometry.bounds.width()
            << "\n"
            << "  height = "
            << internalGeometry.bounds.height()
            << "\n"
            << "  depth  = "
            << internalGeometry.bounds.depth()
            << "\n\n";

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

        //---------------------------------------------------------------------
        // Visualization grid
        //---------------------------------------------------------------------

        const BoundingBox visualizationBox =
            makeVisualizationBox(
                internalGeometry.bounds,
                spacing);


        const std::size_t nx =
            static_cast<std::size_t>(
                std::ceil(
                    visualizationBox.width() /
                    spacing));

        const std::size_t ny =
            static_cast<std::size_t>(
                std::ceil(
                    visualizationBox.height() /
                    spacing));

        const std::size_t nz =
            static_cast<std::size_t>(
                std::ceil(
                    visualizationBox.depth() /
                    spacing));


        const Point origin{
            visualizationBox.min[0] +
                0.5 * spacing,
            visualizationBox.min[1] +
                0.5 * spacing,
            visualizationBox.min[2] +
                0.5 * spacing
        };


        std::cout
            << "\nCell values:\n"
            << "  0.0 = dry\n"
            << "  0.5 = boundary\n"
            << "  1.0 = interior\n\n";


        //---------------------------------------------------------------------
        // Output directories
        //---------------------------------------------------------------------

        const std::filesystem::path outputDirectory =
            std::filesystem::path(
                GEOMETRY_TEST_OUTPUT_DIR) /
            "test_stl_contains_vtks";


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

        const std::vector<double> internalScalar =
            buildScalar(
                internalGeometry,
                visualizationBox,
                spacing);


        writeVTK3DScalar(
            internalDirectory /
                "geometry.vtk",
            nx,
            ny,
            nz,
            origin,
            spacing,
            internalScalar);


        writeVTK3DBoundingBox(
            internalDirectory /
                "bounding_box.vtk",
            internalGeometry.bounds);


        //---------------------------------------------------------------------
        // External
        //---------------------------------------------------------------------

        const std::vector<double> externalScalar =
            buildScalar(
                externalGeometry,
                visualizationBox,
                spacing);


        writeVTK3DScalar(
            externalDirectory /
                "geometry.vtk",
            nx,
            ny,
            nz,
            origin,
            spacing,
            externalScalar);


        writeVTK3DBoundingBox(
            externalDirectory /
                "bounding_box.vtk",
            externalGeometry.bounds);


        writeVTK3DBoundingBox(
            externalDirectory /
                "visualization_box.vtk",
            visualizationBox);


        //---------------------------------------------------------------------
        // Done
        //---------------------------------------------------------------------

        std::cout
            << "Grid dimensions:\n"
            << "  nx = "
            << nx
            << "\n"
            << "  ny = "
            << ny
            << "\n"
            << "  nz = "
            << nz
            << "\n\n";


        std::cout
            << "VTK output written to:\n"
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