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

    if(cellType ==
       CellType::Interior)
    {
        return 1.0;
    }

    throw std::runtime_error(
        "Invalid STL cell type.");
}


BoundingBox makeOpenBox(
    const BoundingBox& bounds)
{
    const double paddingX =
        0.1 * bounds.width();

    const double paddingY =
        0.1 * bounds.height();

    const double paddingZ =
        0.1 * bounds.depth();


    BoundingBox openBox;

    openBox.min =
    {
        bounds.min[0] - paddingX,
        bounds.min[1] - paddingY,
        bounds.min[2] - paddingZ
    };

    openBox.max =
    {
        bounds.max[0] + paddingX,
        bounds.max[1] + paddingY,
        bounds.max[2] + paddingZ
    };


    return openBox;
}


std::vector<double> makeScalar(
    const std::vector<CellType>& cellTypes)
{
    std::vector<double> scalar(
        cellTypes.size());


#pragma omp parallel for schedule(static)
    for(std::ptrdiff_t index = 0;
        index <
            static_cast<std::ptrdiff_t>(
                cellTypes.size());
        ++index)
    {
        const std::size_t cellID =
            static_cast<std::size_t>(
                index);

        scalar[cellID] =
            cellTypeScalar(
                cellTypes[cellID]);
    }


    return scalar;
}


void writeGeometry(
    const STLGeometry& geometry,
    const double spacing,
    const std::filesystem::path& directory)
{
    std::vector<CellType> cellTypes;

    std::vector<double> boundaryX;
    std::vector<double> boundaryY;
    std::vector<double> boundaryZ;


    geometry.interiorAreaAnalysis(
        spacing,
        cellTypes,
        boundaryX,
        boundaryY,
        boundaryZ);


    const BoundingBox& domainBounds =
        geometry.flowType == FlowType::Internal
            ? geometry.bounds
            : geometry.openBox;


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


    const Point origin =
    {
        domainBounds.min[0] + 0.5 * spacing,
        domainBounds.min[1] + 0.5 * spacing,
        domainBounds.min[2] + 0.5 * spacing
    };


    std::cout
        << (geometry.flowType == FlowType::Internal
                ? "Internal"
                : "External")
        << " grid dimensions:\n"
        << "  nx = " << nx << "\n"
        << "  ny = " << ny << "\n"
        << "  nz = " << nz << "\n"
        << "  boundary cells = "
        << boundaryX.size()
        << "\n\n";


    const std::vector<double> scalar =
        makeScalar(
            cellTypes);


    writeVTK3DScalar(
        directory /
            "geometry.vtk",
        nx,
        ny,
        nz,
        origin,
        spacing,
        scalar);


    writeVTK3DBoundingBox(
        directory /
            "bounding_box.vtk",
        geometry.bounds);


    if(geometry.flowType ==
       FlowType::External)
    {
        writeVTK3DBoundingBox(
            directory /
                "open_box.vtk",
            geometry.openBox);
    }
}

} // namespace


int main()
{
    try
    {
        std::cout
            << "========================================\n"
            << "STL Interior Area Analysis Test\n"
            << "========================================\n\n";


        const std::filesystem::path stlFile =
            std::filesystem::path(
                GEOMETRY_TEST_SOURCE_DIR) /
            "stl_tests" /
            "stls" /
            "smooth_irregular_branched_channel.stl";


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
                internalTopology));


        const BoundingBox openBox =
            makeOpenBox(
                internalGeometry.bounds);


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
            FlowType::External,
            &openBox);


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
            << "\n\n"
            << "External open box:\n"
            << "  min    = ("
            << openBox.min[0] << ", "
            << openBox.min[1] << ", "
            << openBox.min[2] << ")\n"
            << "  max    = ("
            << openBox.max[0] << ", "
            << openBox.max[1] << ", "
            << openBox.max[2] << ")\n\n";


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
            "test_stl_vtks";

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


        writeGeometry(
            internalGeometry,
            spacing,
            internalDirectory);

        writeGeometry(
            externalGeometry,
            spacing,
            externalDirectory);


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
