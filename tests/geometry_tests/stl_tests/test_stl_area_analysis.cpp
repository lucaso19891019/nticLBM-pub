#include "geometry_analysis.hpp"
#include "stl_geometry.hpp"
#include "stl_reader.hpp"
#include "stl_validator.hpp"
#include "vtk_output.hpp"
#include "vtk_output_3d.hpp"

#include <cmath>
#include <exception>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <utility>

using namespace ntic::lbm::geometry;
using namespace ntic::lbm::geometry::test;

namespace
{

void printBox(const char* title, const BoundingBox& box)
{
    std::cout
        << title << ":\n"
        << "  min = (" << box.min[0] << ", " << box.min[1] << ", " << box.min[2] << ")\n"
        << "  max = (" << box.max[0] << ", " << box.max[1] << ", " << box.max[2] << ")\n\n";
}

void writeAnalysis(
    const STLGeometry& geometry,
    const double spacing,
    const std::filesystem::path& directory)
{
    GeometryAnalysis3D<ntic::lbm::lattice::LatticeType::D3Q19> analysis(spacing);
    geometry.analysis(analysis);

    std::cout
        << "Analysis domain:\n"
        << "  min = (" << analysis.domain.min[0] << ", "
        << analysis.domain.min[1] << ", " << analysis.domain.min[2] << ")\n"
        << "  max = (" << analysis.domain.max[0] << ", "
        << analysis.domain.max[1] << ", " << analysis.domain.max[2] << ")\n"
        << "Grid dimensions:\n"
        << "  nx = " << analysis.nx << '\n'
        << "  ny = " << analysis.ny << '\n'
        << "  nz = " << analysis.nz << "\n\n";

    writeVTK3DCellScalar(
        directory / "geometry.vtk",
        analysis.nx, analysis.ny, analysis.nz,
        analysis.domain.min, analysis.gridSpacing, analysis.scalar);

    writeVTK3DBoundingBox(
        directory / "bounding_box.vtk", geometry.bounds);
}

} // namespace

int main()
{
    try
    {
        std::cout
            << "========================================\n"
            << "STL Geometry Area Analysis Test\n"
            << "========================================\n\n";

        const std::filesystem::path stlFile =
            std::filesystem::path(GEOMETRY_TEST_SOURCE_DIR) /
            "stl_tests" / "stls" /
            "smooth_irregular_branched_channel.stl";

        std::cout << "STL file:\n  " << stlFile << "\n\n";

        auto internalData = ntic::lbm::stl::read(stlFile.string());
        ntic::lbm::stl::FacetTopology internalTopology;
        ntic::lbm::stl::validate(internalData, "full", internalTopology);
        const STLGeometry internalGeometry(
            std::move(internalTopology), FlowType::Internal);

        // Match the primitive tests: use an explicit open box and show it.
        const double padding = 1.0;
        const BoundingBox openBox{
            {internalGeometry.bounds.min[0] - padding,
             internalGeometry.bounds.min[1] - padding,
             internalGeometry.bounds.min[2] - padding},
            {internalGeometry.bounds.max[0] + padding,
             internalGeometry.bounds.max[1] + padding,
             internalGeometry.bounds.max[2] + padding}
        };

        auto externalData = ntic::lbm::stl::read(stlFile.string());
        ntic::lbm::stl::FacetTopology externalTopology;
        ntic::lbm::stl::validate(externalData, "full", externalTopology);
        const STLGeometry externalGeometry(
            std::move(externalTopology), FlowType::External, &openBox);

        printBox("Bounding box", internalGeometry.bounds);
        printBox("Open box", openBox);

        double spacing = 0.0;
        std::cout << "Enter grid spacing: ";
        std::cin >> spacing;
        if(!std::cin)
            throw std::runtime_error("Failed to read grid spacing.");
        if(!(spacing > 0.0) || !std::isfinite(spacing))
            throw std::runtime_error("Grid spacing must be positive and finite.");

        const std::filesystem::path outputDirectory =
            std::filesystem::path(GEOMETRY_TEST_OUTPUT_DIR) /
            "test_stl_area_analysis_vtks";

        recreateOutputDirectory(outputDirectory);

        const std::filesystem::path internalDirectory =
            outputDirectory / "internal_vtks";
        const std::filesystem::path externalDirectory =
            outputDirectory / "external_vtks";

        std::filesystem::create_directories(internalDirectory);
        std::filesystem::create_directories(externalDirectory);

        std::cout << "\nInternal flow:\n";
        writeAnalysis(internalGeometry, spacing, internalDirectory);

        std::cout << "\nExternal flow:\n";
        writeAnalysis(externalGeometry, spacing, externalDirectory);

        writeVTK3DBoundingBox(
            externalDirectory / "open_box.vtk", openBox);

        std::cout
            << "\nVTK output written to:\n"
            << outputDirectory << "\n";

        return 0;
    }
    catch(const std::exception& exception)
    {
        std::cerr << "Error: " << exception.what() << "\n";
        return 1;
    }
}
