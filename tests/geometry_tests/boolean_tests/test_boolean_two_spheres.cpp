#include "boolean_geometry.hpp"
#include "vtk_output_3d.hpp"

#include <cmath>
#include <exception>
#include <filesystem>
#include <iostream>
#include <stdexcept>

using namespace ntic::lbm::geometry;
using namespace ntic::lbm::geometry::test;

int main()
{
    try
    {
        std::cout << "Boolean two-sphere obstacle: external & external\n\nSphere centers (-1.5,0,0), (1.5,0,0); both radius 2.5.\n\n";
        const BoundingBox openBox{{-7.0, -5.0, -5.0}, {7.0, 5.0, 5.0}};
        const Sphere a({-1.5, 0.0, 0.0}, 2.5, FlowType::External, &openBox);
        const Sphere b({1.5, 0.0, 0.0}, 2.5, FlowType::External, &openBox);
        const auto geometry = a & b;
        double spacing = 0.0;
        std::cout << "Enter grid spacing: ";
        std::cin >> spacing;
        if(!std::cin || !(spacing > 0.0) || !std::isfinite(spacing))
            throw std::runtime_error("Grid spacing must be positive and finite.");

        GeometryAnalysis3D<ntic::lbm::lattice::LatticeType::D3Q19> analysis(spacing);
        geometry.analysis(analysis);

        const std::filesystem::path outputDirectory =
            std::filesystem::path(GEOMETRY_TEST_OUTPUT_DIR) / "test_boolean_two_spheres_vtks";
        recreateOutputDirectory(outputDirectory);

        writeVTK3DCellScalar(
            outputDirectory / "geometry.vtk",
            analysis.nx, analysis.ny, analysis.nz,
            analysis.domain.min, analysis.gridSpacing, analysis.scalar);
        writeVTK3DBoundingBox(outputDirectory / "bounding_box.vtk", geometry.bounds);
        writeVTK3DBoundingBox(outputDirectory / "open_box.vtk", openBox);
        writeVTKBoundaryPoints(outputDirectory / "boundary_q.vtk", analysis);
        std::cout << "\nScalar cells: " << analysis.scalar.size()
                  << "\nBoundary entries: " << analysis.boundaryX.size()
                  << "\nVTK output written to:\n" << outputDirectory << '\n';
        return 0;
    }
    catch(const std::exception& exception)
    {
        std::cerr << "Error: " << exception.what() << '\n';
        return 1;
    }
}