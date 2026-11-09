#pragma once

#include "geometry_analysis.hpp"

#include <algorithm>
#include <array>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

namespace ntic::lbm::geometry::test
{

// Output only Q intersection points. Coincident points with different
// source Face/link records remain separate, preserving provenance.
template<typename Analysis, typename Position, typename Direction>
void writeVTKBoundaryPoints(
    const std::filesystem::path& path,
    const Analysis& analysis,
    Position position,
    Direction direction)
{
    using Lattice = typename Analysis::Lattice;
    constexpr std::size_t links = Lattice::nStencils - 1;
    struct IntersectionPoint
    {
        std::array<double,3> position;
        double q;
        std::size_t face, link;
    };
    if(analysis.q.size() != analysis.boundaryX.size() * links ||
       analysis.boundaryOffsets.size() != analysis.nBoundaries + 1)
        throw std::runtime_error("Incomplete boundary analysis for VTK output.");

    std::vector<IntersectionPoint> points;
    double minimumQ = std::numeric_limits<double>::infinity();
    double maximumQ = 0.0;
    for(std::size_t face = 0; face < analysis.nBoundaries; ++face)
    {
        const auto first = points.size();
        for(std::size_t entry = analysis.boundaryOffsets[face];
            entry < analysis.boundaryOffsets[face + 1]; ++entry)
        {
            const auto start = position(entry);
            for(std::size_t link = 0; link < links; ++link)
            {
                const double q = analysis.q[entry * links + link];
                if(q == -1.0) continue;
                if(!(q > 0.0 && q <= 1.0))
                    throw std::runtime_error("Invalid Q value in boundary output.");
                minimumQ = std::min(minimumQ, q);
                maximumQ = std::max(maximumQ, q);
                auto end = start;
                const auto delta = direction(link + 1);
                for(std::size_t d = 0; d < 3; ++d)
                    end[d] += q * analysis.gridSpacing * delta[d];
                points.push_back({end, q, face, link + 1});
            }
        }
        std::cout << "Face " << face << ": "
                  << analysis.boundaryOffsets[face + 1] - analysis.boundaryOffsets[face]
                  << " entries, " << points.size() - first << " valid Q links\n";
    }

    if(!points.empty())
        std::cout << "Q range: [" << minimumQ << ", " << maximumQ << "]\n";

    std::ofstream output(path);
    if(!output) throw std::runtime_error("Cannot open " + path.string());
    output << std::setprecision(17)
           << "# vtk DataFile Version 3.0\nBoundary Q intersections\nASCII\nDATASET POLYDATA\n"
           << "POINTS " << points.size() << " double\n";
    for(const auto& point : points)
        output << point.position[0] << ' ' << point.position[1] << ' '
               << point.position[2] << '\n';
    output << "VERTICES " << points.size() << ' ' << 2 * points.size() << '\n';
    for(std::size_t i = 0; i < points.size(); ++i)
        output << "1 " << i << '\n';
    output << "POINT_DATA " << points.size() << "\nSCALARS Q double 1\nLOOKUP_TABLE default\n";
    for(const auto& point : points) output << point.q << '\n';
    output << "SCALARS FaceID int 1\nLOOKUP_TABLE default\n";
    for(const auto& point : points) output << point.face << '\n';
    output << "SCALARS LinkID int 1\nLOOKUP_TABLE default\n";
    for(const auto& point : points) output << point.link << '\n';
    if(!output) throw std::runtime_error("Failed to write " + path.string());
}

template<lattice::LatticeType Type>
void writeVTKBoundaryPoints(
    const std::filesystem::path& path, const GeometryAnalysis2D<Type>& analysis)
{
    using Lattice = typename GeometryAnalysis2D<Type>::Lattice;
    writeVTKBoundaryPoints(path, analysis,
        [&](std::size_t entry) { return std::array<double,3>{
            analysis.boundaryX[entry], analysis.boundaryY[entry], 0.0}; },
        [](std::size_t link) { return std::array<int,3>{
            Lattice::ex[link], Lattice::ey[link], 0}; });
}

template<lattice::LatticeType Type>
void writeVTKBoundaryPoints(
    const std::filesystem::path& path, const GeometryAnalysis3D<Type>& analysis)
{
    using Lattice = typename GeometryAnalysis3D<Type>::Lattice;
    writeVTKBoundaryPoints(path, analysis,
        [&](std::size_t entry) { return std::array<double,3>{
            analysis.boundaryX[entry], analysis.boundaryY[entry], analysis.boundaryZ[entry]}; },
        [](std::size_t link) { return std::array<int,3>{
            Lattice::ex[link], Lattice::ey[link], Lattice::ez[link]}; });
}

} // namespace ntic::lbm::geometry::test
