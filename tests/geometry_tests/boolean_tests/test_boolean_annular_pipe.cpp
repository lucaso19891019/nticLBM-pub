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
        std::cout << "Boolean annular pipe: internal & external\n\nCoaxial Z cylinders: outer radius 5 / length 12; inner radius 2 / length 16.\n\n";
        const BoundingBox openBox{{-6.0, -6.0, -9.0}, {6.0, 6.0, 9.0}};
        const Cylinder outer({0.0, 0.0, 0.0}, 5.0, 12.0, Axis::Z, FlowType::Internal);
        const Cylinder inner({0.0, 0.0, 0.0}, 2.0, 16.0, Axis::Z, FlowType::External, &openBox);
        const auto geometry = outer & inner;
        double spacing = 0.0;
        std::cout << "Enter grid spacing: ";
        std::cin >> spacing;
        if(!std::cin || !(spacing > 0.0) || !std::isfinite(spacing))
            throw std::runtime_error("Grid spacing must be positive and finite.");

        GeometryAnalysis3D<ntic::lbm::lattice::LatticeType::D3Q19> analysis(spacing);
        geometry.analysis(analysis);

        const std::filesystem::path outputDirectory =
            std::filesystem::path(GEOMETRY_TEST_OUTPUT_DIR) / "test_boolean_annular_pipe_vtks";
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