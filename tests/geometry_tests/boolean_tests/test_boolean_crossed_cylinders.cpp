#include "boolean_geometry.hpp"
#include "vtk_output_3d.hpp"

#include <cmath>
#include <exception>
#include <filesystem>
#include <iostream>
#include <stdexcept>

using namespace ntic::lbm::geometry;


int main()
{
    try
    {
        std::cout << "Boolean crossed cylinders: internal & internal\n\nCentered X and Z cylinders, both radius 3 / length 10.\n\n";
        const Cylinder a({0.0, 0.0, 0.0}, 3.0, 10.0, Axis::X, FlowType::Internal);
        const Cylinder b({0.0, 0.0, 0.0}, 3.0, 10.0, Axis::Z, FlowType::Internal);
        const auto geometry = a & b;
        double spacing = 0.0;
        std::cout << "Enter grid spacing: ";
        std::cin >> spacing;
        if(!std::cin || !(spacing > 0.0) || !std::isfinite(spacing))
            throw std::runtime_error("Grid spacing must be positive and finite.");

        GeometryAnalysis3D<ntic::lbm::lattice::LatticeType::D3Q19> analysis(spacing);
        geometry.analysis(analysis);

        const std::filesystem::path outputDirectory =
            std::filesystem::path(GEOMETRY_TEST_OUTPUT_DIR) / "test_boolean_crossed_cylinders_vtks";
        recreateOutputDirectory(outputDirectory);

        writeVTK3DCellScalar(
            outputDirectory / "geometry.vtk",
            analysis.nx, analysis.ny, analysis.nz,
            analysis.domain.min, analysis.gridSpacing, analysis.scalar);
        writeVTK3DBoundingBox(outputDirectory / "bounding_box.vtk", geometry.bounds);

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
