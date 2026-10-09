#include "cell_type.hpp"
#include "geometry_analysis.hpp"
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
#include <string>
#include <utility>
#include <vector>

using namespace ntic::lbm::geometry;
using namespace ntic::lbm::geometry::test;

namespace
{
using Analysis = GeometryAnalysis3D<ntic::lbm::lattice::LatticeType::D3Q19>;

void verifyAnalysis(const Analysis& analysis)
{
    const std::size_t count = analysis.nx * analysis.ny * analysis.nz;
    if(count == 0 || analysis.scalar.size() != count)
        throw std::runtime_error("Invalid STL scalar grid dimensions.");
    if(analysis.boundaryX.size() != analysis.boundaryY.size() ||
       analysis.boundaryX.size() != analysis.boundaryZ.size())
        throw std::runtime_error("Boundary coordinate arrays have different sizes.");

    std::size_t boundaryCount = 0;
    for(const double value : analysis.scalar)
    {
        if(value != static_cast<double>(CellType::Dry) &&
           value != static_cast<double>(CellType::Boundary) &&
           value != static_cast<double>(CellType::Interior))
            throw std::runtime_error("Unknown STL cell classification.");
        if(value == static_cast<double>(CellType::Boundary))
            ++boundaryCount;
    }
    if(boundaryCount == 0 || boundaryCount != analysis.boundaryX.size())
        throw std::runtime_error("Boundary count does not match scalar classification.");

    // Every exported boundary center must correspond to a Boundary grid cell.
    const auto coordinateIndex = [&](double x, int axis, std::size_t extent)
    {
        const double normalized = (x - analysis.domain.min[axis]) /
                                  analysis.gridSpacing - 0.5;
        const double nearest = std::round(normalized);
        if(!std::isfinite(normalized) ||
           std::abs(normalized - nearest) > 1.0e-7 ||
           nearest < 0.0 || nearest >= static_cast<double>(extent))
            throw std::runtime_error("Exported boundary coordinate is not a cell center.");
        return static_cast<std::size_t>(nearest);
    };
    std::vector<unsigned char> seen(count, 0);
    for(std::size_t n=0; n<boundaryCount; ++n)
    {
        const std::size_t i=coordinateIndex(analysis.boundaryX[n],0,analysis.nx);
        const std::size_t j=coordinateIndex(analysis.boundaryY[n],1,analysis.ny);
        const std::size_t k=coordinateIndex(analysis.boundaryZ[n],2,analysis.nz);
        const std::size_t id=i+analysis.nx*(j+analysis.ny*k);
        if(seen[id] || analysis.scalar[id] != static_cast<double>(CellType::Boundary))
            throw std::runtime_error("Duplicate or incorrectly classified boundary center.");
        seen[id]=1;
    }
}

void writeAnalysis(const STLGeometry& geometry,
                   double spacing,
                   const std::filesystem::path& directory)
{
    Analysis analysis(spacing);
    geometry.analysis(analysis);
    verifyAnalysis(analysis);

    std::cout << "Grid: " << analysis.nx << " x " << analysis.ny
              << " x " << analysis.nz
              << ", boundary centers: " << analysis.boundaryX.size() << '\n';

    writeVTK3DCellScalar(directory / "geometry.vtk",
                         analysis.nx, analysis.ny, analysis.nz,
                         analysis.domain.min, analysis.gridSpacing, analysis.scalar);
    writeVTK3DBoundingBox(directory / "bounding_box.vtk", geometry.bounds);
}
}

int main(int argc, char** argv)
{
    try
    {
        if(argc > 2)
            throw std::invalid_argument("Usage: test_stl_area_analysis [grid_spacing]");

        const auto stlFile = std::filesystem::path(GEOMETRY_TEST_SOURCE_DIR) /
                             "stl_tests" / "stls" /
                             "smooth_irregular_branched_channel.stl";

        auto internalData=ntic::lbm::stl::read(stlFile.string());
        ntic::lbm::stl::FacetTopology internalTopology;
        ntic::lbm::stl::validate(internalData, "full", internalTopology);
        STLGeometry internalGeometry(std::move(internalTopology), FlowType::Internal);

        const double padding=1.0;
        const BoundingBox openBox{
            {internalGeometry.bounds.min[0]-padding,
             internalGeometry.bounds.min[1]-padding,
             internalGeometry.bounds.min[2]-padding},
            {internalGeometry.bounds.max[0]+padding,
             internalGeometry.bounds.max[1]+padding,
             internalGeometry.bounds.max[2]+padding}
        };

        auto externalData=ntic::lbm::stl::read(stlFile.string());
        ntic::lbm::stl::FacetTopology externalTopology;
        ntic::lbm::stl::validate(externalData, "full", externalTopology);
        STLGeometry externalGeometry(std::move(externalTopology),
                                     FlowType::External, &openBox);

        // Noninteractive default; explicit spacing remains available via argv.
        const double spacing = argc == 2
            ? std::stod(argv[1])
            : internalGeometry.bounds.width()/24.0;
        if(!(spacing > 0.0) || !std::isfinite(spacing))
            throw std::invalid_argument("Grid spacing must be positive and finite.");

        const auto outputDirectory =
            std::filesystem::path(GEOMETRY_TEST_OUTPUT_DIR) / "test_stl_vtks";
        recreateOutputDirectory(outputDirectory);
        const auto internalDirectory=outputDirectory/"internal_vtks";
        const auto externalDirectory=outputDirectory/"external_vtks";
        std::filesystem::create_directories(internalDirectory);
        std::filesystem::create_directories(externalDirectory);

        writeAnalysis(internalGeometry, spacing, internalDirectory);
        writeAnalysis(externalGeometry, spacing, externalDirectory);
        std::cout << "STL area analysis test passed.\n";
        return 0;
    }
    catch(const std::exception& e)
    {
        std::cerr << "STL area analysis test failed: " << e.what() << '\n';
        return 1;
    }
}
