#include "boolean_geometry.hpp"
#include "vtk_output_2d.hpp"

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
        std::cout << "Boolean annulus: internal - internal\n\nConcentric circles: outer radius 5, inner radius 2.\n\n";
        const Circle outer({0.0, 0.0}, 5.0, FlowType::Internal);
        const Circle inner({0.0, 0.0}, 2.0, FlowType::Internal);
        const auto geometry = outer - inner;
        double spacing = 0.0;
        std::cout << "Enter grid spacing: ";
        std::cin >> spacing;
        if(!std::cin || !(spacing > 0.0) || !std::isfinite(spacing))
            throw std::runtime_error("Grid spacing must be positive and finite.");

        GeometryAnalysis2D<ntic::lbm::lattice::LatticeType::D2Q9> analysis(spacing);
        geometry.analysis(analysis);
        if(inner.flowType != FlowType::Internal ||
           geometry.right.flowType != FlowType::Internal)
            throw std::runtime_error("Difference modified the original flow type.");
        const std::filesystem::path outputDirectory =
            std::filesystem::path(GEOMETRY_TEST_OUTPUT_DIR) / "test_boolean_annulus_vtks";
        recreateOutputDirectory(outputDirectory);

        writeVTK2DCellScalar(
            outputDirectory / "geometry.vtk",
            analysis.nx, analysis.ny,
            analysis.domain.min, analysis.gridSpacing, analysis.scalar);
        writeVTK2DBoundingBox(outputDirectory / "bounding_box.vtk", geometry.bounds);

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
