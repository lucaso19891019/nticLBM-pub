#pragma once

#include "bounding_box.hpp"
#include "point.hpp"

#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <vector>


namespace ntic::lbm::geometry::test
{

template <typename FluidFunction>
void writeVTK3DGeometry(
    const std::filesystem::path& path,
    const BoundingBox& domainBounds,
    const double spacing,
    FluidFunction isFluid)
{
    if(spacing <= 0.0)
    {
        throw std::invalid_argument(
            "Grid spacing must be positive.");
    }


    const std::size_t nx =
        static_cast<std::size_t>(
            std::floor(
                domainBounds.width() /
                spacing)) +
        1;

    const std::size_t ny =
        static_cast<std::size_t>(
            std::floor(
                domainBounds.height() /
                spacing)) +
        1;

    const std::size_t nz =
        static_cast<std::size_t>(
            std::floor(
                domainBounds.depth() /
                spacing)) +
        1;


    const std::size_t pointCount =
        nx * ny * nz;


    std::vector<int> fluid(
        pointCount,
        0);


    const std::size_t xySize =
        nx * ny;


#pragma omp parallel for schedule(static)
    for(std::ptrdiff_t index = 0;
        index <
            static_cast<std::ptrdiff_t>(
                pointCount);
        ++index)
    {
        const std::size_t pointIndex =
            static_cast<std::size_t>(
                index);


        const std::size_t k =
            pointIndex /
            xySize;

        const std::size_t remainder =
            pointIndex %
            xySize;

        const std::size_t j =
            remainder /
            nx;

        const std::size_t i =
            remainder %
            nx;


        const double x =
            domainBounds.min[0] +
            static_cast<double>(i) *
            spacing;

        const double y =
            domainBounds.min[1] +
            static_cast<double>(j) *
            spacing;

        const double z =
            domainBounds.min[2] +
            static_cast<double>(k) *
            spacing;


        const Point point{
            x,
            y,
            z
        };


        fluid[pointIndex] =
            isFluid(point)
                ? 1
                : 0;
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
        << "nticLBM 3D geometry test\n"
        << "ASCII\n"
        << "DATASET STRUCTURED_POINTS\n"
        << "DIMENSIONS "
        << nx << " "
        << ny << " "
        << nz << "\n"
        << "ORIGIN "
        << domainBounds.min[0] << " "
        << domainBounds.min[1] << " "
        << domainBounds.min[2] << "\n"
        << "SPACING "
        << spacing << " "
        << spacing << " "
        << spacing << "\n"
        << "POINT_DATA "
        << pointCount << "\n"
        << "SCALARS fluid int 1\n"
        << "LOOKUP_TABLE default\n";


    for(const int value :
        fluid)
    {
        output
            << value
            << "\n";
    }
}


//=============================================================================
// 3D cell scalar field
//=============================================================================

inline void writeVTK3DCellScalar(
    const std::filesystem::path& path,
    const std::size_t nx,
    const std::size_t ny,
    const std::size_t nz,
    const Point& origin,
    const double spacing,
    const std::vector<double>& scalar)
{
    if(spacing <= 0.0)
    {
        throw std::invalid_argument(
            "Grid spacing must be positive.");
    }

    const std::size_t cellCount =
        nx * ny * nz;

    if(scalar.size() !=
       cellCount)
    {
        throw std::invalid_argument(
            "3D scalar size does not match "
            "the VTK cell dimensions.");
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
        << "nticLBM 3D cell scalar field\n"
        << "ASCII\n"
        << "DATASET STRUCTURED_POINTS\n"
        << "DIMENSIONS "
        << nx + 1 << " "
        << ny + 1 << " "
        << nz + 1 << "\n"
        << "ORIGIN "
        << origin[0] << " "
        << origin[1] << " "
        << origin[2] << "\n"
        << "SPACING "
        << spacing << " "
        << spacing << " "
        << spacing << "\n"
        << "CELL_DATA "
        << cellCount << "\n"
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


inline void writeVTK3DBoundingBox(
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


    const double xmin =
        box.min[0];

    const double ymin =
        box.min[1];

    const double zmin =
        box.min[2];

    const double xmax =
        box.max[0];

    const double ymax =
        box.max[1];

    const double zmax =
        box.max[2];


    output
        << "# vtk DataFile Version 3.0\n"
        << "nticLBM 3D bounding box\n"
        << "ASCII\n"
        << "DATASET POLYDATA\n"
        << "POINTS 8 double\n"

        << xmin << " "
        << ymin << " "
        << zmin << "\n"

        << xmax << " "
        << ymin << " "
        << zmin << "\n"

        << xmax << " "
        << ymax << " "
        << zmin << "\n"

        << xmin << " "
        << ymax << " "
        << zmin << "\n"

        << xmin << " "
        << ymin << " "
        << zmax << "\n"

        << xmax << " "
        << ymin << " "
        << zmax << "\n"

        << xmax << " "
        << ymax << " "
        << zmax << "\n"

        << xmin << " "
        << ymax << " "
        << zmax << "\n"

        << "LINES 12 36\n"

        << "2 0 1\n"
        << "2 1 2\n"
        << "2 2 3\n"
        << "2 3 0\n"

        << "2 4 5\n"
        << "2 5 6\n"
        << "2 6 7\n"
        << "2 7 4\n"

        << "2 0 4\n"
        << "2 1 5\n"
        << "2 2 6\n"
        << "2 3 7\n";
}

} // namespace ntic::lbm::geometry::test