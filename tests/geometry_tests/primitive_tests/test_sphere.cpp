#include "sphere.hpp"
#include "geometry_analysis.hpp"
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

void writeAnalysis(
    const Sphere& geometry,
    const double spacing,
    const std::filesystem::path& directory)
{
    GeometryAnalysis3D<ntic::lbm::lattice::LatticeType::D3Q19> analysis(
        spacing);

    geometry.analysis(
        analysis);

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
        geometry.boundingBox());
}

} // namespace


int main()
{
    try
    {
        std::cout
            << "========================================\n"
            << "Sphere Geometry Test\n"
            << "========================================\n\n";

        const BoundingBox openBox{
            {-1.0, -2.0, 0.0},
            {9.0, 8.0, 10.0}
        };

        const Sphere internalSphere{
            {4.0, 3.0, 5.0},
            2.0,
            FlowType::Internal
        };

        const Sphere externalSphere{
            {4.0, 3.0, 5.0},
            2.0,
            FlowType::External,
            &openBox
        };

        const BoundingBox boundingBox =
            internalSphere.boundingBox();

        std::cout
            << "Sphere:\n"
            << "  center = ("
            << internalSphere.center[0] << ", "
            << internalSphere.center[1] << ", "
            << internalSphere.center[2] << ")\n"
            << "  radius = "
            << internalSphere.radius
            << "\n\n";

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
            "test_sphere_vtks";

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

        writeAnalysis(
            internalSphere,
            spacing,
            internalDirectory);

        writeAnalysis(
            externalSphere,
            spacing,
            externalDirectory);

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
