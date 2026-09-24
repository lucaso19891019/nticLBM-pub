#pragma once

#include "bounding_box.hpp"
#include "point.hpp"

#include <cstddef>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <vector>


namespace ntic::lbm::geometry::test
{

//=============================================================================
// 2D scalar field
//=============================================================================

inline void writeVTK2DScalar(
    const std::filesystem::path& path,
    const std::size_t nx,
    const std::size_t ny,
    const Point& origin,
    const double spacing,
    const std::vector<double>& scalar)
{
    if(spacing <= 0.0)
    {
        throw std::invalid_argument(
            "Grid spacing must be positive.");
    }


    const std::size_t pointCount =
        nx * ny;


    if(scalar.size() !=
       pointCount)
    {
        throw std::invalid_argument(
            "2D scalar size does not match "
            "the VTK grid dimensions.");
    }


    std::ofstream output(path);

    if(!output)
    {
        throw std::runtime_error(
            "Failed to open VTK output file: " +
            path.string());
    }


    output
        << "# vtk DataFile Version 3.0\n"
        << "nticLBM 2D scalar field\n"
        << "ASCII\n"
        << "DATASET STRUCTURED_POINTS\n"
        << "DIMENSIONS "
        << nx << " "
        << ny << " "
        << 1 << "\n"
        << "ORIGIN "
        << origin[0] << " "
        << origin[1] << " "
        << origin[2] << "\n"
        << "SPACING "
        << spacing << " "
        << spacing << " "
        << 1.0 << "\n"
        << "POINT_DATA "
        << pointCount << "\n"
        << "SCALARS scalar double 1\n"
        << "LOOKUP_TABLE default\n";


    for(const double value :
        scalar)
    {
        output
            << value
            << "\n";
    }
}


//=============================================================================
// 2D bounding box
//=============================================================================

inline void writeVTK2DBoundingBox(
    const std::filesystem::path& path,
    const BoundingBox& box)
{
    std::ofstream output(path);

    if(!output)
    {
        throw std::runtime_error(
            "Failed to open VTK output file: " +
            path.string());
    }


    output
        << "# vtk DataFile Version 3.0\n"
        << "nticLBM 2D bounding box\n"
        << "ASCII\n"
        << "DATASET POLYDATA\n"
        << "POINTS 4 double\n"

        << box.min[0] << " "
        << box.min[1] << " "
        << box.min[2] << "\n"

        << box.max[0] << " "
        << box.min[1] << " "
        << box.min[2] << "\n"

        << box.max[0] << " "
        << box.max[1] << " "
        << box.min[2] << "\n"

        << box.min[0] << " "
        << box.max[1] << " "
        << box.min[2] << "\n"

        << "LINES 4 12\n"
        << "2 0 1\n"
        << "2 1 2\n"
        << "2 2 3\n"
        << "2 3 0\n";
}

} // namespace ntic::lbm::geometry::test