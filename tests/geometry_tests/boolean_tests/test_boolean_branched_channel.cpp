#include "boolean_geometry.hpp"
#include "stl_reader.hpp"
#include "stl_validator.hpp"
#include "vtk_output_3d.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <exception>
#include <filesystem>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

using namespace ntic::lbm::geometry;
using namespace ntic::lbm::geometry::test;

namespace
{

std::vector<Cylinder> makeGuides(const STLGeometry& channel)
{
    if(channel.boundaryFeatures.nPlanarFaces != 7)
        throw std::runtime_error("Expected seven planar ports in the branched-channel STL.");

    std::vector<Cylinder> guides;
    for(const auto& face : channel.boundaryFeatures.planarFaces)
    {
        Point low{std::numeric_limits<double>::infinity(),
                  std::numeric_limits<double>::infinity(),
                  std::numeric_limits<double>::infinity()};
        Point high{-low[0], -low[1], -low[2]};
        std::vector<std::size_t> vertices;
        Point normal{};
        for(std::size_t f = 0; f < channel.boundaryFeatures.facetFaceIDs.size(); ++f)
        {
            if(channel.boundaryFeatures.facetFaceIDs[f] != face.faceID) continue;
            const auto& facet = channel.topology.facetGeometry[f];
            for(std::size_t d = 0; d < 3; ++d) normal[d] += facet.area * facet.normal[d];
            for(const auto vertex : channel.topology.geometry.facetVertexIDs[f])
            {
                vertices.push_back(vertex);
                const auto& p = channel.topology.geometry.vertices[vertex];
                for(std::size_t d = 0; d < 3; ++d)
                {
                    low[d] = std::min(low[d], p[d]);
                    high[d] = std::max(high[d], p[d]);
                }
            }
        }
        if(vertices.empty()) throw std::runtime_error("Port has no facets.");
        std::size_t axis = 0;
        for(std::size_t d = 1; d < 3; ++d)
            if(std::abs(normal[d]) > std::abs(normal[axis])) axis = d;
        if(!(std::abs(normal[axis]) > 0.0))
            throw std::runtime_error("Port has no oriented normal.");
        for(std::size_t d = 0; d < 3; ++d)
            if(d != axis && std::abs(normal[d]) > 1e-8 * std::abs(normal[axis]))
                throw std::runtime_error("Port normal is not aligned with a Cylinder axis.");

        Point center{};
        for(std::size_t d = 0; d < 3; ++d) center[d] = (low[d] + high[d]) / 2.0;
        double radius = 0.0;
        for(const auto vertex : vertices)
        {
            const auto& p = channel.topology.geometry.vertices[vertex];
            double squared = 0.0;
            for(std::size_t d = 0; d < 3; ++d)
                if(d != axis) squared += (p[d] - center[d]) * (p[d] - center[d]);
            radius = std::max(radius, std::sqrt(squared));
        }
        // Cover the irregular port contour and overlap the original channel.
        radius *= 1.02;
        const double length = 4.0 * radius;
        const double overlap = radius;
        // STLGeometry normals point toward the solid. At this internal
        // channel's port cap, that is the outward extension direction.
        center[axis] += (normal[axis] > 0.0 ? 1.0 : -1.0) * (length / 2.0 - overlap);
        const std::array<Axis,3> axes{Axis::X, Axis::Y, Axis::Z};
        guides.emplace_back(center, radius, length, axes[axis], FlowType::Internal);
        std::cout << "Port Face " << face.faceID << " -> guide " << guides.size()
                  << ": axis " << axis << ", center (" << center[0] << ", "
                  << center[1] << ", " << center[2] << "), radius " << radius
                  << ", length " << length << ", overlap " << overlap << '\n';
    }
    if(guides.size() != 7) throw std::runtime_error("Expected seven guide cylinders.");
    return guides;
}

} // namespace

int main()
{
    try
    {
        std::cout << "Boolean branched channel: nested internal unions\n\n";
        const auto stlFile = std::filesystem::path(GEOMETRY_TEST_SOURCE_DIR) /
            "stl_tests" / "stls" / "smooth_irregular_branched_channel.stl";
        std::cout << "STL file: " << stlFile << '\n';
        auto data = ntic::lbm::stl::read(stlFile.string());
        ntic::lbm::stl::FacetTopology topology;
        ntic::lbm::stl::validate(data, "full", topology);
        STLGeometry channel(std::move(topology), FlowType::Internal);
        auto guides = makeGuides(channel);

        double spacing = 0.0;
        std::cout << "Enter grid spacing: ";
        std::cin >> spacing;
        if(!std::cin || !(spacing > 0.0) || !std::isfinite(spacing))
            throw std::runtime_error("Grid spacing must be positive and finite.");

        auto step1 = std::move(channel) | guides[0];
        auto step2 = std::move(step1) | guides[1];
        auto step3 = std::move(step2) | guides[2];
        auto step4 = std::move(step3) | guides[3];
        auto step5 = std::move(step4) | guides[4];
        auto step6 = std::move(step5) | guides[5];
        const auto geometry = std::move(step6) | guides[6];

        GeometryAnalysis3D<ntic::lbm::lattice::LatticeType::D3Q19> analysis(spacing);
        geometry.analysis(analysis);

        const auto outputDirectory = std::filesystem::path(GEOMETRY_TEST_OUTPUT_DIR) /
            "test_boolean_branched_channel_vtks";
        recreateOutputDirectory(outputDirectory);
        writeVTK3DCellScalar(outputDirectory / "geometry.vtk",
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
