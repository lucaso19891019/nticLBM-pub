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
void writeVTK2DGeometry(
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


    const std::size_t pointCount =
        nx * ny;


    std::vector<int> fluid(
        pointCount,
        0);


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


        const std::size_t i =
            pointIndex % nx;

        const std::size_t j =
            pointIndex / nx;


        const double x =
            domainBounds.min[0] +
            static_cast<double>(i) *
            spacing;

        const double y =
            domainBounds.min[1] +
            static_cast<double>(j) *
            spacing;


        const Point point{
            x,
            y,
            0.0
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
        << "nticLBM 2D geometry test\n"
        << "ASCII\n"
        << "DATASET STRUCTURED_POINTS\n"
        << "DIMENSIONS "
        << nx << " "
        << ny << " 1\n"
        << "ORIGIN "
        << domainBounds.min[0] << " "
        << domainBounds.min[1] << " 0\n"
        << "SPACING "
        << spacing << " "
        << spacing << " 1\n"
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
        << box.min[1] << " 0\n"

        << box.max[0] << " "
        << box.min[1] << " 0\n"

        << box.max[0] << " "
        << box.max[1] << " 0\n"

        << box.min[0] << " "
        << box.max[1] << " 0\n"

        << "LINES 4 12\n"

        << "2 0 1\n"
        << "2 1 2\n"
        << "2 2 3\n"
        << "2 3 0\n";
}

} // namespace ntic::lbm::geometry::test