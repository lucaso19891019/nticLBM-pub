#include "geometry_analysis.hpp"
#include "stl_geometry.hpp"
#include "stl_reader.hpp"
#include "stl_validator.hpp"
#include "vtk_output.hpp"
#include "vtk_output_3d.hpp"

#include <exception>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <utility>


using namespace ntic::lbm::geometry;
using namespace ntic::lbm::geometry::test;


namespace
{

void printAnalysis(
    const GeometryAnalysis3D<void>& analysis)
{
    std::cout
        << "Analysis domain:\n"
        << "  min = ("
        << analysis.domain.min[0] << ", "
        << analysis.domain.min[1] << ", "
        << analysis.domain.min[2] << ")\n"
        << "  max = ("
        << analysis.domain.max[0] << ", "
        << analysis.domain.max[1] << ", "
        << analysis.domain.max[2] << ")\n"
        << "Grid dimensions:\n"
        << "  nx = "
        << analysis.nx << '\n'
        << "  ny = "
        << analysis.ny << '\n'
        << "  nz = "
        << analysis.nz << "\n\n";
}


void writeAnalysis(
    const STLGeometry& geometry,
    const double spacing,
    const std::filesystem::path& directory)
{
    GeometryAnalysis3D<void> analysis(
        spacing);

    geometry.analysis(
        analysis);

    printAnalysis(
        analysis);

    writeVTK3DCellScalar(
        directory /
            "geometry.vtk",
        analysis.nx,
        analysis.ny,
        analysis.nz,
        analysis.domain.min,
        analysis.gridSpacing,
        analysis.scalar);

    writeVTK3DBoundingBox(
        directory /
            "bounding_box.vtk",
        geometry.bounds);
}

} // namespace


int main()
{
    try
    {
        std::cout
            << "========================================\n"
            << "STL Area Analysis Test\n"
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
                internalTopology),
            FlowType::Internal);

        const double padding =
            1.0;

        const BoundingBox openBox{
            {
                internalGeometry.bounds.min[0] - padding,
                internalGeometry.bounds.min[1] - padding,
                internalGeometry.bounds.min[2] - padding
            },
            {
                internalGeometry.bounds.max[0] + padding,
                internalGeometry.bounds.max[1] + padding,
                internalGeometry.bounds.max[2] + padding
            }
        };

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
            << externalGeometry.openBox.min[0] << ", "
            << externalGeometry.openBox.min[1] << ", "
            << externalGeometry.openBox.min[2] << ")\n"
            << "  max    = ("
            << externalGeometry.openBox.max[0] << ", "
            << externalGeometry.openBox.max[1] << ", "
            << externalGeometry.openBox.max[2] << ")\n"
            << "  width  = "
            << externalGeometry.openBox.width()
            << "\n"
            << "  height = "
            << externalGeometry.openBox.height()
            << "\n"
            << "  depth  = "
            << externalGeometry.openBox.depth()
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

        std::cout
            << "\nInternal analysis:\n";

        writeAnalysis(
            internalGeometry,
            spacing,
            internalDirectory);

        std::cout
            << "External analysis:\n";

        writeAnalysis(
            externalGeometry,
            spacing,
            externalDirectory);

        std::cout
            << "VTK output directory:\n  "
            << outputDirectory
            << "\n\n";

        std::cout
            << "STL area analysis test passed.\n";

        return 0;
    }
    catch(const std::exception& error)
    {
        std::cerr
            << "STL area analysis test failed: "
            << error.what()
            << '\n';

        return 1;
    }
}
